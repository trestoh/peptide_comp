// pep_target_comp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>

int main(int argc, char* argv[])
{
	std::ifstream perc_in(argv[1]);
	std::ifstream rts_in(argv[2]);

	std::string output = argv[3] + std::string("/match_results.txt");

	std::ofstream perc_out;
	perc_out.open(output, std::ios_base::app);

	std::string junk;
	std::string line;
	std::string peptide;
	double q_val;
	int scan;
	
	std::vector<std::string> perc_peptides;
	std::vector<std::string> rts_peptides;

	std::getline(perc_in, line);

	while (std::getline(perc_in, line))
	{
		std::istringstream iss(line);
		if (!(iss >> junk >> scan >> junk >> junk >> junk >> junk >> junk >> q_val >> junk >> junk >> peptide)) { break; } // error

		if (q_val < .01)
			perc_peptides.push_back(peptide);

		//perc_out << peptide << std::endl;
	}

	std::string xcorr1, xcorr2, xcorr3, deltacn;
	std::getline(rts_in, line);
	std::istringstream isss(line);
	isss >> xcorr1 >> xcorr2 >> xcorr3 >> deltacn;

	while (std::getline(rts_in, line))
	{
		std::istringstream iss(line);
		if (!(iss >> peptide)) { break; }

		rts_peptides.push_back(peptide);

	}

	std::unordered_map<std::string, bool> in_rts;
	std::unordered_map<std::string, bool> in_perc;
	
	for (int i = 0; i < rts_peptides.size(); i++)
	{
		in_rts.insert({ rts_peptides[i], true });
	}

	for (int i = 0; i < perc_peptides.size(); i++)
	{
		in_perc.insert({ perc_peptides[i], true });
	}

	int wrong_rts = 0;
	int missed_targets = 0;

	for (int i = 0; i < rts_peptides.size(); i++)
	{
		std::string find_me = rts_peptides[i];
		auto search = in_perc.find(find_me);
		if (search == in_perc.end())
		{
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
			//perc_out << "Peptide: " << perc_peptides[i] << " in Percolator but not RTS" << std::endl;
			missed_targets++;
		}
	}

	perc_out << "Parmeters: (" << xcorr1 << ", " << xcorr2 << ", " << xcorr3 << ", " << deltacn << ") with "
		<< wrong_rts << " false positives over Percolator and " << missed_targets << " targets missed from percolator" << std::endl;

	//std::sort(peptides.begin(), peptides.end());

	//for (int i = 0; i < peptides.size(); i++)
	//{
	//	perc_out << peptides[i] << std::endl;
	//}

	return 0;
}
