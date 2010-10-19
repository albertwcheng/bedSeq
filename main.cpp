/***************************************************************************
 Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include "StringUtil.h"
#include "RandomAccessFile.h"
#include <vector>
#include "Nucleic.h"
#include <deque>
using namespace std;

#define  BEDLINE_BUFFER_SIZE 1024*1024

#define BED_FORMAT 1
#define EBED_FORMAT 2

#define FORWARD '+'
#define REVERSE '-'

int printUsageAndExit(const char*programName)
{
	cerr<<"Usage: "<<programName<<" <seqDir> <bedfilename> <bedfiletype:bed|ebed> [options] > <ofilename>"<<endl;
	cerr<<"Description: \tappend sequence to bed files (bed format for the feature coordinate, ebed format for feature block coordinates)"<<endl;
	cerr<<"\t\tseqDir contains per chromosome <chr>.seq which is one line raw sequence without header. convert fa to seq file using convertFaToPureSeq.py"<<endl;
	cerr<<"Options:"<<endl;
	cerr<<"--output-fasta\toutput fasta file instead of appending seq to bed file"<<endl;
	cerr<<"--use-coord-as-name\tuse coordinate as name of sequence instead of the name field of bed file. Only valid when --output-fasta option is specified"<<endl;
	cerr<<"--use-block-coord-as-name\tuse block coordinate as name. Only when --output-fasta is specified and format is ebed"<<endl;
	cerr<<"--print-OK\tBy default, only error messages are printed to stderr. specifying --print-OK ask to print OK to stderr even when succesful to allow for efficient handling of piping from a master program"<<endl;
	return EXIT_FAILURE;

}

inline string getSeqFromRaf(RandomAccessFile* raf,int start0, int end1,char strand)
{
	string forward=raf->get(start0,end1);
	if(strand==FORWARD)
		return forward;
	
	return reverse_complement(forward);
}

int main(int argc, const char**argv){
	
	char *bedLineBuffer =new char[BEDLINE_BUFFER_SIZE];
	
	bool outputFasta=false;
	bool useCoordAsName=false;
	bool useBlockCoordAsName=false;
	bool printOK=false;
	
	const char* programName=argv[0];
	
	
	if(argc<4){
		return printUsageAndExit(programName);
	}
	
	for(int ai=4;ai<argc;ai++)
	{
		if(!strcmp(argv[ai],"--output-fasta"))
			outputFasta=true;
		else if(!strcmp(argv[ai],"--use-coord-as-name"))
			useCoordAsName=true;
		else if(!strcmp(argv[ai],"--use-block-coord-as-name"))
			useBlockCoordAsName=true;
		else if(!strcmp(argv[ai],"--print-OK"))
			printOK=true;
		else {
			cerr<<"unknown option "<<argv[ai]<<endl;
			return printUsageAndExit(programName);
		}

	}
	
	string seqDir=argv[1];
	string bedfilename=argv[2];
	string bedfiletype=argv[3];
	
	int bedformat=0;
	
	map<string,RandomAccessFile*> rafStreams;
	
	if(bedfiletype=="ebed")
		bedformat=EBED_FORMAT;
	else if (bedfiletype=="bed")
		bedformat=BED_FORMAT;
	else
	{
		cerr<<"Unknown bed format:"<<bedfiletype<<".abort"<<endl;
		return printUsageAndExit(programName);
	}
	
	
	if(useBlockCoordAsName && useCoordAsName ){
		cerr<<"both --use-block-coord-as-name and --use-coord-as-name are specified. These two options are mutually exclusive"<<endl;
		return printUsageAndExit(programName);
	}
	
	if( (useBlockCoordAsName || useCoordAsName ) && ( !outputFasta) ){
		cerr<<"--use-block-coord-as-name and --use-coord-as-name are valid only when --output-fasta is specified"<<endl;
		return printUsageAndExit(programName);
	}
	
	if( useBlockCoordAsName && bedformat!=EBED_FORMAT){
		cerr<<"--use-block-coord-as-name is only valid when format is ebed"<<endl;
		return printUsageAndExit(programName);
	}
	
	
	//now open and read bed file
	ifstream bedfil(bedfilename.c_str());
	
	strcpy(bedLineBuffer,"\0");
	vector<string> fields;
	vector<string> blockSizes;
	vector<string> blockStarts;
	
	int lino=0;
	
	string thisChr="";
	RandomAccessFile* thisRaf=NULL;
	
	while(bedfil.good())
	{
		if(!bedfil.good())
			break;
		
		lino++;
		bedfil.getline(bedLineBuffer,BEDLINE_BUFFER_SIZE);
		StringUtil::split(bedLineBuffer,"\t",fields);
		
		if(strlen(bedLineBuffer)==0)
			continue;
		
		if(fields.size()<3)
		{
			cerr<<"Ignored: invalid number of fields for line "<<lino<<" ["<<bedLineBuffer<<"]"<<endl;
			continue;
		}
		
		string chrom(fields[0]);
	
		if (chrom!=thisChr)
		{
			//switch raf!!
			//
			//
			map<string,RandomAccessFile*>::iterator rafI=rafStreams.find(chrom);
			if( rafI==rafStreams.end())
			{
				thisRaf=new RandomAccessFile(seqDir+"/"+chrom+".seq");
				//register
				rafStreams.insert(map<string,RandomAccessFile*>::value_type(chrom,thisRaf));
			}
			else {
				thisRaf=rafI->second;
			}

		}
		
		int featureStart0=StringUtil::atoi(fields[1]);
		int featureEnd1=StringUtil::atoi(fields[2]);
		
		char strand=FORWARD;
		
		if (fields.size()>=6) //get strand info
		{
			
			if(fields[5]=="-")
				strand=REVERSE;
			else 
				strand=FORWARD;
			

		}

		
		if(bedformat==BED_FORMAT)
		{
			//simple!
			string seq=getSeqFromRaf(thisRaf,featureStart0,featureEnd1,strand);
			if(outputFasta)
			{
				string name;
				if(useCoordAsName || fields.size()<4)
				{	name=chrom+":"+StringUtil::str(featureStart0+1)+"-"+StringUtil::str(featureEnd1);
					
				}
				else {
					name=fields[3];
				}
				
				cout<<">"<<name;
				cout<<endl;
				cout<<seq;
				cout<<endl;
				if(printOK)
					cerr<<"OK"<<endl;

			}
			else {
				cout<<bedLineBuffer;
				cout<<"\t";
				cout<<seq;
				cout<<endl;
				if(printOK)
					cerr<<"OK"<<endl;
			}

			
		}
		else {
			//this is ebed format
			if(fields.size()<12)
			{
				cerr<<"Ignored: invalid number of fields for line "<<lino<<" "<<bedLineBuffer<<endl;
				continue;
			}
			
			int blockCount=StringUtil::atoi(fields[9]);
			
			string seq="";
			StringUtil::split(fields[10],",",blockSizes);
			StringUtil::split(fields[11],",",blockStarts);
			
			if (blockSizes.size()<blockCount || blockStarts.size()<blockCount)
			{
				cerr<<"Ignored: invalid number of blockStarts and blockEnds given the blockCount "<<lino<<" "<<bedLineBuffer<<endl;
				continue;
			}

			for (int blockI=0;blockI<blockCount;blockI++)
			{
				int blockRelStart=StringUtil::atoi(blockStarts[blockI]);
				int blockSize=StringUtil::atoi(blockSizes[blockI]);
				int getStart0=featureStart0+blockRelStart;
				int getEnd1=getStart0+blockSize;
				seq+=getSeqFromRaf(thisRaf,getStart0,getEnd1,FORWARD); //always get forward first then convert at end
			}
			
			
			if(outputFasta){
				string name;
				if(useCoordAsName)
				{	name=chrom+":"+StringUtil::str(featureStart0+1)+"-"+StringUtil::str(featureEnd1);
					
				}
				else if(useBlockCoordAsName){
					
					
					deque<string> coordblocks;
					
					for (int blockI=0;blockI<blockCount;blockI++)
					{
						int blockRelStart=StringUtil::atoi(blockStarts[blockI]);
						int blockSize=StringUtil::atoi(blockSizes[blockI]);
						int getStart0=featureStart0+blockRelStart;
						int getEnd1=getStart0+blockSize;
						string coordstring=chrom+":"+StringUtil::str(getStart0+1)+"-"+StringUtil::str(getEnd1);
						if(strand==FORWARD)
							coordblocks.push_back(coordstring);
						else {
							coordblocks.push_front(coordstring);
						}

					}
					
					name=StringUtil::join< deque<string>, deque<string>::const_iterator >(coordblocks,",");
				}
				else {
					name=fields[3];
				}
				
				cout<<">"<<name;
				cout<<endl;
				cout<<seq;
				cout<<endl;
				if(printOK)
					cerr<<"OK"<<endl;
				
			}else{
				
			
				cout<<bedLineBuffer;
				cout<<"\t";
				if(strand==FORWARD)
					cout<<seq;
				else {
					cout<<reverse_complement(seq);
				}
				cout<<endl;
				if(printOK)
					cerr<<"OK"<<endl;
			}
			
			
			
		}

		
	}
	
	bedfil.close();
	
	//close all rafs
	for(map<string,RandomAccessFile*>::iterator i=rafStreams.begin();i!=rafStreams.end();i++)
	{
		i->second->close();
		delete i->second;
	}
	
	delete[] bedLineBuffer;
	
	return EXIT_SUCCESS;
}