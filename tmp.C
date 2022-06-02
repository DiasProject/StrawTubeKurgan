void tmp(){

    	ifstream fin("result.txt"); // открыли файл для чтения
  	ofstream fout("cppstudio.txt"); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
    	
	TFile hfile("Event.root","RECREATE");
 
	Int_t counter, id, zerocounter;
	counter=0;   

   	TTree *tree = new TTree("T","tree");
	tree->Branch("counter",&counter,"counter/I");

	bool ready=0;


	while(!fin.eof( )) {


	fin >> id;
	if(id < 9999){
	counter++;
	cout<<counter<<endl;
	}
	else {
	//	if(counter>0){
	//	cout<< counter<<endl;
	//	fout<< counter<<endl;
		tree->Fill();
		counter=0;
		
		}
	}
	bool tmp=0;
	counter=0;	
	while(!fin.eof( )) {

	
	fin >> id;
	if(id > 9999){
		ready=1;
		if(ready==tmp)counter++;
		else tmp=ready;
	}
	else {
		ready = 0;
		tmp = 0;
	}
	//if(ready && id > 9999) counter++;
	
	}
	cout<<counter<<endl;
  	tree->Write();
    	fin.close(); // закрываем файл
	fout.close();
	hfile.Close();

}
