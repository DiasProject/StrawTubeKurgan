#include "KurganDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4SDManager.hh"

#include "comments.hh"
#include "math.h"


using namespace std; 

KurganDetectorConstruction::KurganDetectorConstruction()
: G4VUserDetectorConstruction(),
  fSensetiveDetector(nullptr),
  fSensetiveDetector2(nullptr),
  fScoringVolume(nullptr)
{ }

KurganDetectorConstruction::~KurganDetectorConstruction()
{ 
}

G4VPhysicalVolume* KurganDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = false; //true;

  //////////////////////////////*** WORLD ***//////////////////////////////
  
  //G4double world_sizeXY = 750.*cm;	//1.2*env_sizeXY;
  //G4double world_sizeZ  = 750.*cm;	//1.2*env_sizeZ;
  G4double world_sizeXY = 200.*m;	//1.2*env_sizeXY;
  G4double world_sizeZ  = 200.*m;	//1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       1.0*world_sizeXY, 1.0*world_sizeXY, 1.0*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //////////////////////////////*** VACUUM ***//////////////////////////////
                     
  G4double atomicNumber = 1.;   
  G4double massOfMole = 1.008*g/mole;   
  G4double density = 1.e-25*g/cm3;   
  G4double temperature = 300*kelvin;   
  G4double pressure = 10.E-8*bar; 
  G4Material* Vacuum =      
	new G4Material("interGalactic", 
			atomicNumber,
         massOfMole, 
			density, 
			kStateGas,
			temperature, 
			pressure);

  //////////////////////////////*** ENVELOPE ***//////////////////////////////
  
  //G4double env_sizeXY = 600. * cm; 
  //G4double env_sizeZ = 600. * cm;
  G4double env_sizeXY = 200.*m; 
  G4double env_sizeZ = 200.*m;
  G4Material* env_mat = Vacuum;

  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        1.0*env_sizeXY, 1.0*env_sizeXY, 1.0*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
      
  //////////////////////////////*** GAS ***//////////////////////////////


   			/////***Elements***/////
   G4double density_CO2 = 0.001913*g/cm3; 
   G4Material* CO2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE"); // density = 0.00184212(g/cm3)
   
   G4double density_Ar = 0.001784*g/cm3; 
   G4Material* Ar = nist->FindOrBuildMaterial("G4_Ar"); // density = 0.00166201(g/cm3)
       		
			/////***Gas in straw tubes***/////  
   // 30% CO2 + 70% Ar
   G4double density_gas = 0.3 * density_CO2 + 0.7 * density_Ar;//0.001848*g/cm3
   G4Material* Straw_Gas = new G4Material("Straw_Gas_1", density_gas, 2);   
   Straw_Gas->AddMaterial(CO2, 30*perCent);   
   Straw_Gas->AddMaterial(Ar, 70*perCent);   
   //cout<<"density_(0.3*C02+0.7*Ar)="<<density_gas/(g/cm3)<<"(g/cm3)"<<endl; 
   
   //////////////////////////////*** STRAW CHAMBER ***//////////////////////////////
   //материалы для строу трубки
   G4Material* straw_mat_MY = nist->FindOrBuildMaterial("G4_MYLAR"); //density(g/cm^3) = 1.4,l(eV) = 166  
   G4Material* straw_mat_Cu = nist->FindOrBuildMaterial("G4_Cu"); //density(g/cm^3) = 8.96,  l(eV) = 322, thikness = 50 nm
   //G4double density_Au = 19.31*g/cm3;
   G4Material* straw_mat_Au = nist->FindOrBuildMaterial("G4_Au"); //density(g/cm^3) = 19.31, l(eV) = 790, thikness = 20 nm
   //материалы для Анода
   //G4double density_W = 19.25*g/cm3;
   G4Material* anod_mat_W = nist->FindOrBuildMaterial("G4_W"); //density(g/cm^3) = 19.25, Diameter = 20 mkm
    
	//толщина слоев и покрытий строу трубки    
   G4double maylar_film = 0.02*mm;//36.0*pow(10.0,-3)*mm;
   G4double Cu_coating = 50.0*nm;//pow(10.0,-7)*cm;
   G4double Au_coating = 20.0*nm;//pow(10.0,-7)*cm;
   

  //химический состав почвы 
/*  G4double density_O = 0.00133151*g/cm3;
  G4double density_Si = 2.33*g/cm3;
  G4double density_Al = 2.699*g/cm3;
  G4double density_Fe = 7.874*g/cm3;
  G4double density_C = 2*g/cm3;
  G4double density_Ca = 1.55*g/cm3;
  G4double density_K = 0.862*g/cm3;
  G4double density_Na = 0.971*g/cm3;
  G4double density_Mg = 1.74*g/cm3;
  G4double density_N = 0.0011652*g/cm3;
*/

  G4Material* O = nist->FindOrBuildMaterial("G4_O");
  G4Material* Si = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Al = nist->FindOrBuildMaterial("G4_Al");
  G4Material* Fe = nist->FindOrBuildMaterial("G4_Fe");
  G4Material* C = nist->FindOrBuildMaterial("G4_C");
  G4Material* Ca = nist->FindOrBuildMaterial("G4_Ca");
  G4Material* K = nist->FindOrBuildMaterial("G4_K");
  G4Material* Na = nist->FindOrBuildMaterial("G4_Na");
  G4Material* Mg = nist->FindOrBuildMaterial("G4_Mg"); 
  G4Material* N = nist->FindOrBuildMaterial("G4_N"); 
  
//  G4double density_Pb = 11.35*g/cm3;
  G4Material* Pb_mat = nist->FindOrBuildMaterial("G4_Pb"); 

//  G4double density_O2 = 0.00142897*g/cm3;
//  G4Material* O2 = new G4Material("O2", density_O2, 1);
//  O2->AddElement(O2, 2);

  G4double density_Soil = 5.51*g/cm3;
  G4Material* Soil_mat = new G4Material("Soil", density_Soil, 10);
  Soil_mat->AddMaterial(O, 100*perCent);
  Soil_mat->AddMaterial(Si, 33*perCent);
  Soil_mat->AddMaterial(Al, 8.11*perCent);
  Soil_mat->AddMaterial(Fe, 3.8*perCent);
  Soil_mat->AddMaterial(C, 2*perCent);
  Soil_mat->AddMaterial(Ca, 1.37*perCent);
  Soil_mat->AddMaterial(K, 1.36*perCent);
  Soil_mat->AddMaterial(Na, 0.63*perCent);
  Soil_mat->AddMaterial(Mg, 0.63*perCent);
  Soil_mat->AddMaterial(N, 0.1*perCent);

/*  G4double density_P = 2.2*g/cm3;

  G4double density_OxideKalium = 2.32*g/cm3;
  G4double density_OxidePhosfor = 2.39*g/cm3;
  G4Material* K = nist->FindOrBuildMaterial("G4_K");
  G4Material* O = nist->FindOrBuildMaterial("G4_O");
  G4Material* P = nist->FindOrBuildMaterial("G4_P") 
  G4Material* N = nist->FindOrBuildMaterial("G4_N") 


  G4Material* K2O = new G4Material("OxideKalium", density_OxideKalium, 2); 
  K2O->AddElement(K, 2);
  K2O->AddElement(O, 1);
  
  G4Material* P2O5 = new G4Material("OxidePhosfor", density_OxidePhosfor, 2);
  P2O5->AddElement(P, 2);
  P2O5->AddElement(O, 5);

  G4double density_Soil = 5.51*g/cm3;
  G4Material* Soil = new G4Material("Soil", density_Soil, 3);
  Soil->AddMaterial(K2O, 30*perCent);
  Soil->AddMaterial(P2O5, 30*perCent);
  Soil->AddMaterial(N, 30*perCent);
*/

   G4double high_Soil = 100.*cm;
   G4double width_Soil = 100.*cm;
   G4double length_Soil = 35.*m;
   G4Box* Kurgan_Box = new G4Box("Kurgan", width_Soil, high_Soil, length_Soil);
    
    //fSolidVolume
   G4LogicalVolume* Kurgan_Volume =  new G4LogicalVolume(Kurgan_Box, /*Straw_Gas*/Pb_mat /*straw_mat_MY*/, "Kurgan");
    
   G4ThreeVector position_Kurgan = G4ThreeVector(0.*cm, 0.*cm, length_Soil+5.2*m);
    new G4PVPlacement(0,                       //no rotation
                    position_Kurgan,         //at (0,0,0)
                    Kurgan_Volume,                //its logical volume
                    "Kurgan_Box",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


	//геометрические характеристики Анода
   G4double anod_Au = 10.0*nm;
   G4double D_anod = 0.02*mm;
   G4double R_anod = D_anod/2;
   G4double R_anod_W = R_anod - anod_Au;
   
	//геометрические характеристики строу трубки
   G4double high = 100.*cm;
   G4double Straw_outer_diameter = 20.0*mm; 
   G4double Straw_outer_R = Straw_outer_diameter/2;
   G4double Straw_R_1 = Straw_outer_R - maylar_film;
   G4double Straw_R_2 = Straw_outer_R - maylar_film - Cu_coating;
   G4double Straw_R_3 = Straw_outer_R - maylar_film - Cu_coating - Au_coating;
   G4double Straw_Dz = high/2;
   G4double Straw_SPhi = 0.*deg;
   G4double Straw_DPhi = 360.*deg;

   G4Tubs* Straw_Maylar = new G4Tubs("Straw", Straw_R_1, Straw_outer_R, Straw_Dz, Straw_SPhi, Straw_DPhi);//слой Майлара
   G4Tubs* Straw_Cu = new G4Tubs("Straw", Straw_R_2, Straw_R_1, Straw_Dz, Straw_SPhi, Straw_DPhi);//покрытие из Меди
   G4Tubs* Straw_Au = new G4Tubs("Straw", Straw_R_3, Straw_R_2, Straw_Dz, Straw_SPhi, Straw_DPhi);//покрытие из Золота
   G4Tubs* Straw_Inner_Gas = new G4Tubs("Straw", R_anod, Straw_R_3-0.1*mm, Straw_Dz, Straw_SPhi, Straw_DPhi);//Газ в строу трубке
   G4Tubs* Straw_Anod_Au = new G4Tubs("Straw", R_anod_W, R_anod, Straw_Dz, Straw_SPhi, Straw_DPhi);//Газ в строу трубке
   G4Tubs* Straw_Anod_W = new G4Tubs("Straw", 0.*mm, R_anod_W, Straw_Dz, Straw_SPhi, Straw_DPhi);//Газ в строу трубке  
   //G4Box* Aperture = new G4Box("Aperture", Straw_Dz+20.*cm, Straw_Dz+20.*cm, 0.001*cm);
   G4Box* Aperture = new G4Box("Aperture", width_Soil, high_Soil, 0.001*cm);
	
#ifdef SHOW_COMMENTS_DetectorConstruction_cc  
	cout<<"Straw_outer_D="<<Straw_outer_diameter/mm<<"mm"<<endl;
	cout<<"Straw_outer_R="<<Straw_outer_R/mm<<"mm"<<endl;
	cout<<"Straw_R_1="<<Straw_R_1/mm<<"mm"<<endl;
	cout<<"Straw_R_2="<<Straw_R_2/mm<<"mm"<<endl;
	cout<<"Straw_R_3="<<Straw_R_3/mm<<"mm"<<endl;
	cout<<"Straw_R_3-0.1="<<(Straw_R_3-0.1*mm)/mm<<"mm"<<endl;
	cout<<"R_anod="<<R_anod/mm<<"mm"<<endl;
	cout<<"Anod_AU="<<anod_Au/mm<<"mm"<<endl;
	cout<<"Anod_W="<<R_anod_W/mm<<"mm"<<endl;
	cout<<"D_effect="<<(2*(Straw_R_3-0.1*mm))/mm<<" mm(current)"<<endl;
	cout<<"D_effect="<<(2*Straw_R_3)/mm<<" mm(before)"<<endl;
#endif

  //первый слой строу трубки из Майлара
   G4LogicalVolume* logicShape1 =                    
     new G4LogicalVolume(Straw_Maylar,         //its solid
                        straw_mat_MY,          //its material
                        "Straw");           //its name
    //второй слой(покрытие 50нм) строу трубки из Меди  
   G4LogicalVolume* logicShape2 =                    
     new G4LogicalVolume(Straw_Cu,         //its solid
                        straw_mat_Cu,          //its material
                        "Straw");           //its name
    //третий слой(покрытие 20нм) строу трубки из Золота   
   G4LogicalVolume* logicShape3 =                    
     new G4LogicalVolume(Straw_Au,         //its solid
                        straw_mat_Au,          //its material
                        "Straw");           //its name
    //Газ внутри строу трубки   
    //fSensetiveVolume
   G4LogicalVolume* logicShape4 =                    
     new G4LogicalVolume(Straw_Inner_Gas,         //its solid
                        Straw_Gas,          //its material
                        "StrawGas");           //its name
    //Анод Au  
   G4LogicalVolume* logicShape5 =                    
     new G4LogicalVolume(Straw_Anod_Au,         //its solid
                        straw_mat_Au,          //its material
                        "Straw");           //its name
    //Анод W   
   G4LogicalVolume* logicShape6 =                    
     new G4LogicalVolume(Straw_Anod_W,         //its solid
                        anod_mat_W,          //its material
                        "Straw");           //its name

    //fSensetiveVolume
   G4LogicalVolume* logicShape7 =                    
     new G4LogicalVolume(Aperture,         //its solid
                        Straw_Gas,          //its material
                        "Aperture");           //its name

    //fSensetiveVolume
   G4LogicalVolume* logicShape8 =                    
     new G4LogicalVolume(Aperture,         //its solid
                        Straw_Gas,          //its material
                        "Aperture");           //its name
    //////////////////////////////*** View X ***//////////////////////////////
    //первая камера, состоит из 4 слоев(XXYY), строу трубок 200 трубки на слой 
    //вторая камера, состоит из 4 слоев(XXYY), строу трубок 200 трубки на слой 
    //третяя камера, состоит из 4 слоев(XXYY), строу трубок 200 трубки на слой 
 
   G4double position_XY = 0.0;
   G4double layer_position_Z = 0.0;
   G4double phi;
   G4int N_layer = 4;																     //количество слоев в камере (X,X,Y,Y)
   G4double d_between_layer = Straw_outer_diameter;			              //расстояние между слоями в камере
   G4double d_between_chamber = 500.0*mm;                            	  //расстояние между камерами
   G4double d_between_straw  = Straw_outer_diameter;				           //расстояние между центрами строу трубок
   //G4double XY_offset = Straw_outer_R;								              //расстояние от пучковой линий до первой строу второго блока
   G4double ChamberPosition_Z = 0.0*mm;
   G4int N_straw_tot = 600;														     //общее количесвто строу трубок
   G4int N_chamber = 3;                                                   //общее количество камер 3
   //G4int Views = 2;															           //количество плоскостей 2 (X,Y)
   G4int N_straws_per_chamber = N_straw_tot/N_chamber;					     //количество строу трубок на каждой камере 
   G4int N_straws_per_layer = N_straws_per_chamber/N_layer;				     //количество строу трубок на каждом слое 
	 G4double offset_layer = Straw_outer_R;                                 //сдвег между слоями по оси XY

#ifdef SHOW_COMMENTS_DetectorConstruction_cc  
	cout<<"N_straw_tot="<<N_straw_tot<<endl;
	//cout<<"Views="<<Views<<endl;
	cout<<"N_straws_per_chamber="<<N_straws_per_chamber<<endl;
	cout<<"N_layer="<<N_layer<<endl;
	cout<<"N_straws_per_layer="<<N_straws_per_layer<<endl;
#endif

/*
   for(G4int chamber=0; chamber<N_chamber; chamber++){
    cout<<"ChamberPosition_Z = "<<ChamberPosition_Z<<endl;    
      for(G4int layer=0; layer<N_layer; layer++){ 
         layer_position_Z = ChamberPosition_Z + d_between_layer * layer;
   
   G4ThreeVector position;
   G4RotationMatrix rot = G4RotationMatrix();
   G4RotationMatrix rot_1 = G4RotationMatrix();
   G4Transform3D transform = G4Transform3D();
         if(layer<N_layer/2){
            
            phi=0.0*deg; 
            rot.rotateX(phi); 
            rot.rotateY(phi); 
            rot.rotateZ(phi); 
            rot_1.rotateX(90.0*deg);

         }else{
            
            phi=90.0*deg; 
            rot.rotateX(0.0*deg); 
            rot.rotateY(phi); 
            rot.rotateZ(0.0*deg); 
            rot_1.rotateY(90.0*deg);

        }
      

      for(G4int i=0; i<N_straws_per_layer; i++){
	   		if(i<25){
					position_XY = i * d_between_straw;	
            if(layer==1)position_XY = position_XY + offset_layer;
            if(layer==3)position_XY = position_XY + offset_layer;
		   	}else if(i>=25 && i<50){
					position_XY = (i-N_straws_per_layer)*d_between_straw;
            if(layer==1)position_XY = position_XY + offset_layer;
            if(layer==3)position_XY = position_XY + offset_layer;
			  }

  #ifdef SHOW_COMMENTS_DetectorConstruction_cc_test 
    G4cout<<i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1<<" "<<position_XY*cos(phi)<<" "<<position_XY*sin(phi)<<" "<<layer_position_Z<<G4endl;
    G4cout<<rot<<" "<<rot_1<<G4endl;
  #endif
    
      position = G4ThreeVector(position_XY*cos(phi), position_XY*sin(phi), layer_position_Z);
		  transform = G4Transform3D(rot_1, position);
	
			new G4PVPlacement(transform,			     //rot, position
                    logicShape1,             //its logical volume
		    			      "Straw",                 //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1,                        //copy number
                    checkOverlaps);          //overlaps checking
			 
      new G4PVPlacement(transform, 			     //rot, position
                    logicShape2,             //its logical volume
		        			  "Straw",                 //its name
					      	  logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1,                       //copy number
                    checkOverlaps);          //overlaps checking

      new G4PVPlacement(transform,			     //rot, position
                    logicShape3,             //its logical volume
			      		    "Straw",                 //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1,              //copy number
                    checkOverlaps);          //overlaps checking

			new G4PVPlacement(transform,			     //rot, position
                    logicShape4,             //its logical volume
					          "Straw",                 //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1,                       //copy number
                    checkOverlaps);          //overlaps checking
           
	   	new G4PVPlacement(transform,			     //rot, position
                    logicShape5,             //its logical volume
					          "Anod",                  //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1,                       //copy number
                    checkOverlaps);          //overlaps checking
		
	   	new G4PVPlacement(transform,			     //rot, position
                    logicShape6,             //its logical volume
					          "Anod",                  //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    i+N_straws_per_layer*layer+chamber*N_layer*N_straws_per_layer+1,                       //copy number
                    checkOverlaps);          //overlaps checking
		}
	}
  ChamberPosition_Z = ChamberPosition_Z + d_between_chamber;   
}
*/
/////////////////sensetive detector////////////////  
   G4ThreeVector position_Aperture = G4ThreeVector(0.*cm, 0.*cm, 0.2*m);
		  new G4PVPlacement(0,			 				    //rotm
			      		position_Aperture,	 			  //position
                logicShape7,           		  //its logical volume
	        		  "Aperture",               	//its name
                logicEnv,                		//its mother  volume
                false,                   		//no boolean operation
                N_straw_tot+1,             	//copy number
                checkOverlaps);          		//overlaps checking

   G4ThreeVector position_Aperture_2 = G4ThreeVector(0.*cm, 0.*cm, 80.2*m);
		  new G4PVPlacement(0,			 				    //rotm
			      		position_Aperture_2,	 			  //position
                logicShape8,           		  //its logical volume
	        		  "Aperture",               	//its name
                logicEnv,                		//its mother  volume
                false,                   		//no boolean operation
                602,             	//copy number
                checkOverlaps);          		//overlaps checking

  fScoringVolume = logicShape4;
  fSensetiveDetector = logicShape7;
  fSensetiveDetector2 = logicShape8;
 
  if(checkOverlaps){cout<<"checkOverlaps!"<<endl;}   
  return physWorld;
}
