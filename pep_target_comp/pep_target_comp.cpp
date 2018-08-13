// pep_target_comp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <iostream>

int main(int argc, char* argv[])
{
	std::ifstream perc_in(argv[1]);
	std::ifstream rts_in(argv[2]);

	std::string output = argv[3] + std::string("/match_summary.txt");


	std::string master_out_name = std::string("C:/dev/workspace/master_match.txt");

	std::ofstream perc_out;
	std::ofstream master_out;

	perc_out.open(output, std::ios_base::app);
	master_out.open(master_out_name, std::ios_base::app);

	if (argc > 3)
	{
		std::string file_name = argv[4];
		std::cout << "File name is: " << file_name << std::endl;
		master_out << file_name << "\t";
	}

	std::string junk;
	std::string line;
	std::string peptide;
	double q_val;
	int scan;
	
	std::vector<std::string> perc_peptides;
	std::vector<std::string> rts_peptides;
	std::vector<int> perc_scan;
	std::vector<int> rts_scan;

	std::getline(perc_in, line);

	while (std::getline(perc_in, line))
	{
		std::istringstream iss(line);
		if (!(iss >> junk >> scan >> junk >> junk >> junk >> junk >> junk >> q_val >> junk >> junk >> peptide)) { break; } // error

		if (q_val < .01)
		{
			perc_peptides.push_back(peptide);
			perc_scan.push_back(scan);
		}
		//perc_out << peptide << std::endl;
	}

	std::string xcorr1, xcorr2, xcorr3, deltacn;
	std::getline(rts_in, line);
	std::istringstream isss(line);
	isss >> xcorr1 >> xcorr2 >> xcorr3 >> deltacn;

	while (std::getline(rts_in, line))
	{
		std::istringstream iss(line);
		if (!(iss >> scan >> peptide)) { break; }

		rts_peptides.push_back(peptide);
		rts_scan.push_back(scan);

	}

	std::string smart_params;
	for (int i = 0; i < xcorr1.length(); i++)
		if (xcorr1[i] != '.')
			smart_params.push_back(xcorr1[i]);

	for (int i = 0; i < xcorr2.length(); i++)
		if (xcorr2[i] != '.')
			smart_params.push_back(xcorr2[i]);

	for (int i = 0; i < xcorr3.length(); i++)
		if (xcorr3[i] != '.')
			smart_params.push_back(xcorr3[i]);

	for (int i = 0; i < deltacn.length(); i++)
		if (deltacn[i] != '.')
			smart_params.push_back(deltacn[i]);

	//perc_out << "Param combo is: " << smart_params << std::endl;

	std::string params = xcorr1 + xcorr2 + xcorr3 + deltacn;
	std::size_t slug = std::hash<std::string>{}(params);

	slug = slug >> 48;

	//std::string tag = std::to_string(slug);

	std::string tag = smart_params;

	std::string false_hit_out = argv[3] + std::string("/false_hits_") + tag + std::string(".txt");
	std::string miss_target_out = argv[3] + std::string("/missed_targets_") + tag + std::string(".txt");

	std::ofstream false_out;
	false_out.open(false_hit_out, std::ios_base::app);

	std::ofstream target_out;
	target_out.open(miss_target_out, std::ios_base::app);

	//perc_out << "Slug is: " << slug << " with byte size: " << sizeof(size_t) << std::endl;

	false_out << "Params are: " << xcorr1 << " " << xcorr2 << " " << xcorr3 << " " << deltacn << std::endl;
	target_out << "Params are: " << xcorr1 << " " << xcorr2 << " " << xcorr3 << " " << deltacn << std::endl;


	std::unordered_map<std::string, int> in_rts;
	std::unordered_map<std::string, int> in_perc;
	
	for (int i = 0; i < rts_peptides.size(); i++)
	{
		in_rts.insert({ rts_peptides[i], rts_scan[i] });
	}

	for (int i = 0; i < perc_peptides.size(); i++)
	{
		in_perc.insert({ perc_peptides[i], perc_scan[i] });
	}

	int wrong_rts = 0;
	int missed_targets = 0;

	for (int i = 0; i < rts_peptides.size(); i++)
	{
		std::string find_me = rts_peptides[i];
		auto search = in_perc.find(find_me);
		if (search == in_perc.end())
		{
			false_out << rts_scan[i] << '\t' << rts_peptides[i] << std::endl;
			//perc_out << "Peptide: " << rts_peptides[i] << " in RTS but not Percolator" << std::endl;
			wrong_rts++;
		}
	}

	for (int i = 0; i < perc_peptides.size(); i++)
	{
		std::string find_me = perc_peptides[i];
		auto search = in_rts.find(find_me);
		if (search == in_rts.end())
		{
			target_out << perc_scan[i] << '\t' << perc_peptides[i] << std::endl;
			//perc_out << "Peptide: " << perc_peptides[i] << " in Percolator but not RTS" << std::endl;
			missed_targets++;
		}
	}

	//perc_out << "Parmeters: (" << xcorr1 << ", " << xcorr2 << ", " << xcorr3 << ", " << deltacn << ") with "
	//	<< wrong_rts << " false positives over Percolator and " << missed_targets << " targets missed from percolator" << std::endl;

	double pct_over = (double)wrong_rts / (double)perc_peptides.size();
	double pct_under = 1 - (double)(perc_peptides.size() - missed_targets) / (double)perc_peptides.size();

	perc_out << xcorr1 << "\t" << xcorr2 << "\t" << xcorr3 << "\t" << deltacn << "\t" << missed_targets 
		<< "\t" << wrong_rts << "\t" << pct_under << "\t" << pct_over << std::endl ;

	master_out << xcorr1 << "\t" << xcorr2 << "\t" << xcorr3 << "\t" << deltacn << "\t" << missed_targets
		<< "\t" << wrong_rts << "\t" << pct_under << "\t" << pct_over << std::endl;


	//std::sort(peptides.begin(), peptides.end());

	//for (int i = 0; i < peptides.size(); i++)
	//{
	//	perc_out << peptides[i] << std::endl;
	//}

	return 0;
}
