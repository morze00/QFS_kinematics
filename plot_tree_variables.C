void plot_tree_variables(TString rootfilename)
{
	TFile* file = TFile::Open(rootfilename);

    if (!file || file->IsZombie()) {

        std::cout << "Cannot open file\n";
        return;
    }
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cout << "Tree not found\n";
        return;
    }

	auto c1 = new TCanvas("c1", "c1", 1200, 900);
	c1->Divide(4, 3);
	
	c1->cd(1);
	tree->Draw("theta_1:theta_2>>h1(180,0,180,180,0,180)","","colz");
		
	c1->cd(2);
	tree->Draw("theta_1>>h2(180,0,180)","","");

	c1->cd(3);
	tree->Draw("theta_2>>h3(180,0,180)","","");
	
	c1->cd(4);
	tree->Draw("th1_cm:th2_cm>>h4(180,0,180,180,0,180)","","colz");
		
	c1->cd(5);
	tree->Draw("th1_cm>>h5(180,0,180)","","");

	c1->cd(6);
	tree->Draw("th2_cm>>h6(180,0,180)","","");
	
	c1->cd(7);
	tree->Draw("Opang>>h7(180,0,180)","","");

	c1->cd(8);
	tree->Draw("Dif_phi>>h8(360,0,360)","","");

	c1->cd(9);
	tree->Draw("phi_1:phi_2>>h9(360,-180,180,360,-180,180)","","colz");
	
	c1->cd(10);
	tree->Draw("E1:E2","","colz");

	c1->cd(11);
	tree->Draw("PBx>>h_PBx_PBy_PBz_rf(200,-500,500)","","");
	TH1F* h_PBx = (TH1F*)gDirectory->Get("h_PBx_PBy_PBz_rf");
	h_PBx->SetLineColor(kRed);

	tree->Draw("PBy>>h_PBy(200,-500,500)", "", "same");
	TH1F *h_PBy = (TH1F *)gDirectory->Get("h_PBy");
	h_PBy->SetLineColor(kBlue);

	tree->Draw("PBz_rf>>h_PBz_rf(200,-500,500)", "", "same");
	TH1F *h_PBz_rf = (TH1F *)gDirectory->Get("h_PBz_rf");
	h_PBz_rf->SetLineColor(kGreen);

	c1->cd(12);
	tree->Draw("theta_B>>h_theta_B(1000,0,180)", "", "");
}
