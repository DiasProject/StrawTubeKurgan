void count()
{
	ifstream input_file;
	ofstream output_file;
	input_file.open("counter_mu.txt");
	output_file.open("output.txt");

	TFile hfile("Event.root","RECREATE");
	TTree *tree = new TTree("T","tree");
	Int_t tmp, counter=0;
	tree->Branch("counter",&counter,"counter/I");

	while(!input_file.eof()){
	input_file>>tmp;
		if(tmp!=1000){
			counter++;
		}else{
		tree->Fill();
		if(counter==0){output_file<<counter<<endl;}
		counter=0;}
	}
	
	tree->Write();
	hfile.Close();
	input_file.close();
	output_file.close();
}
