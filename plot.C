void plot()
{

	gROOT->SetStyle("Plain");
	TTree *tree = new TTree("tree","treelibrated tree");
	Int_t sum = 0;
	Int_t counter;
	tree->Branch("event",&sum, "sum/I");
	// fill some events with random numbers
	ifstream in_file("test.txt");
	while(!in_file.eof())
	{ 
		in_file>>counter;
			sum = counter;
			//cout<<sum<<endl;	
			tree->Fill();
	}

   // now draw some tree variables
	TCanvas *c1 = new TCanvas();
	c1->cd(1);
	tree->Draw("sum");  //energy of det a
	tree->Write();
	in_file.close();
}
