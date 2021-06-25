#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h> 
#include "zlib.h"




#include "BamCommonLibrary.h"

//#include "Function.h"




void 	alignment_List(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	alignment_Coverage(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length);
void 	alignment_ListAndCoverage(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name);
void	strcat_triple(char	*main, char *first, char *second, char *third, int length);

void	usage(){
	printf("Program: bam-amplicon (for BAM file)\n");
	printf("Version: 1.0\n");
	printf("Usage:   bam-amplicon -m <Mode> -b <bam_filename> [Options]\n");
	printf("\n");
	//printf("Mode:	depthdist	about depth's distribution\n");
	printf("Mode:	depthdist	show depths\n");
	printf("	quality		show distribution of quality scores for a specified position\n");
	//printf("	ampsummary	about amplicon's depth and cover rate\n");
	printf("\n");
	printf("	pattern		show frequency of patterns within a specified region\n");
	printf("	del		show frequency of deletions with position info\n");
	printf("	ins		show read sequences with marked insertions\n");
	printf("\n");
	printf("	stat		show BAM statistics\n");
	printf("	length		show distribution of mapped read lenghts\n");
	printf("	ampsummary	show depths and coverage rates within specified amplicons\n");
//	printf("	trim		trim the BAM file using amplicons provided by the BED file\n");
	printf("	trim		show the nearest amplicon for trimming,by using amplicons provided by the BED file\n");
	printf("\n");

	printf("Options:\n");
	printf("	-b, --bam    [FILE]	input BAM file\n");
	printf("	-r, --bed    [FILE]	the range of specified regions from BED file\n");
	printf("	-c, --chr    [STR]	chromosome name\n");
	printf("	-s, --start  [INT]	start position (1-based, i.e. first base is 1)\n");
	printf("	-e, --end    [INT]	end position (1-based, i.e. first base is 1)\n");
//	printf("	-g, --fasta  [FILE]	fasta file (reference file)\n");
	printf("	-t, --target [FILE]	the ragne of amplicon regions from BED file\n");
//	printf("	-f, --filter [TYPE]	filter criteria\n");
//	printf("				[readqual, mapqual]\n");
	printf("	-l, --compact		show only average depth and coverage rate (depthdist mode) \n");
	printf("	-a, --origin		show original sequence (include insertion)\n");
//	printf("	-d, --duplicate		show duplicate\n");
//	printf("	-v, --verbose		show the processed status \n");
//	printf("	-n, --column_name	remove the first row (header)\n");
	printf("	-u, --threshold	[INT]	coverage above the threshold\n");
//	printf("	/* For Ion Torrent Data */\n");
//	printf("	-w, --flow		modify the sequence by flow signals (FZ)\n");
//	printf("	-z, --zm		change the FZ values to ZM values\n");

	printf("	-p, --pattern	[STR]	show the percentage of ALT pattern specified (pattern mode only)\n");
	printf("	-h, --help		show the manual page\n");
	printf("	\n");
}


int main(int argc, char *argv[]){


//New Variable
//
	FILE	*file_bam_i = NULL;
	FILE	*file_bai_i = NULL;
	FILE	*file_bed_i = NULL;
	FILE	*file_fasta_i = NULL;
	FILE	*file_fai_i = NULL;

	char	*mode;
	char	*bam;
	char	*bai_1;
	char	*bai_2;
	char	*fasta;
	char	*fai;
	char	*bed;
	char	*target;

	int	flag_bam;
	int	flag_fasta = 0;
//END

	toolsFlags	ToolsFlags;



	const char* const short_options = "m:b:r:c:s:e:g:t:o:f:p:u:vndhlaMzw";
	struct option long_options[] = {
		{"mode"	,	1,	NULL,	'm'},
		{"bam"	,	1,	NULL,	'b'},
		{"bed"	,	1,	NULL,	'r'},
		{"chr"	,	1,	NULL,	'c'},
		{"start",	1,	NULL,	's'},
		{"end",		1,	NULL,	'e'},
		{"target",	1,	NULL,	't'},
		{"fasta",	1,	NULL,	'g'},
		{"output",	1,	NULL,	'o'},
		{"filter",	1,	NULL,	'f'},
		{"pattern",	1,	NULL,	'p'},
		{"threshold",	1,	NULL,	'u'},
		{"verbose",	0,	NULL,	'v'},
		{"column_name",	0,	NULL,	'n'},
		{"compact",	0,	NULL,	'l'},
		{"origin",	0,	NULL,	'a'},
		{"duplicate",	0,	NULL,	'd'},
		{"mapq",	0,	NULL,	'M'},
		{"zm",		0,	NULL,	'z'},
		{"flow",	0,	NULL,	'w'},
		{"help",	0,	NULL,	'h'}
	};

	int	c;

	if (argc == 1){
		usage();
		return -1;
	}

	flag_bam	= 0;


	//ALL Flag is Zero initially
	memset(&ToolsFlags,0,sizeof(toolsFlags));


	while ((c = getopt_long (argc, argv, short_options, long_options, NULL)) != -1){
		switch (c){
			case 'm':	
				mode = strdup(optarg);	
				//printf("mode:\t%s\n",optarg);	
				break;
			case 'b':
				//Check Bam & Bai
				bam = strdup(optarg);
				bai_1 = calloc(strlen(bam)+6,sizeof(char));
				strcpy(bai_1, optarg);	strcat( bai_1, ".bai");
				bai_2 = strdup(optarg);	bai_2[strlen(bai_2)-1] = 'i';

				if ( (access(bam, R_OK) != -1) && ((access(bai_1, R_OK) != -1) || (access(bai_2, R_OK) != -1))){	
					flag_bam = 1;
				}else {
					if (access(bam, R_OK) == -1){
						printf("bam:\t%s\n",bam);	
						printf("Please Check Bam File.\n");
					}else {
						printf("bai_1:\t%s\n",bai_1);	
						printf("bai_2:\t%s\n",bai_2);
						printf("Please Check Bai File.\n");
					}
					return -1;
				}
				break;
			case 'g':
				fasta	= strdup(optarg);
				fai	= calloc(strlen(fasta)+6,sizeof(char));
				strcpy(fai, optarg);	strcat(fai, ".fai");
				if ( (access(fasta, R_OK) != -1) && (access(fai, R_OK) != -1) ){	
					flag_fasta = 1;
				}else {
					if (access(fasta, R_OK) == -1){
						printf("fasta:\t%s\n", fasta);	
						printf("Please Check Fasta File.\n");
					}else {
						printf("fai:\t%s\n", fai);	
						printf("Please Check Fai File.\n");
					}
				}
				break;
			case 'r':	
				//Check Bed File
				bed = strdup(optarg);
				if ( access(bed, R_OK) != -1){
					//printf("bed:\t%s\n",bed);
					ToolsFlags.flag_bed = 1;
					/*Pattern*/
					ToolsFlags.flag_pattern = 1;
				}else {
					printf("Please Check Bed File.\n");
					return -1;
				}
				break;
			case 'c':
				ToolsFlags.chromosome = strdup(optarg);
				ToolsFlags.flag_chromosome = 1;
				break;
			case 's':	
				ToolsFlags.start = atoi(optarg);
				ToolsFlags.flag_start = 1;
				if (ToolsFlags.start < 1){
					printf("\n[Error!]\tStart Point:[%d] is smaller than 1.\n\n", ToolsFlags.start);
					usage();
					return -1;
				}
				break;
			case 'e':
				ToolsFlags.end	= atoi(optarg);
				ToolsFlags.flag_end = 1;
				//printf("end:\t%s\n",optarg);	
				break;
			case 't':	
				target = strdup(optarg);
				if ( access(target, R_OK) != -1){
					ToolsFlags.flag_target = 1;
				}else {
					printf("Please Check Bed File.\n");
					return -1;
				}
				break;
			case 'o':	printf("output:\t%s\n",optarg);	break;
			case 'f':	
				ToolsFlags.flag_filter = 1;
				//printf("filter:\t%s\n",optarg);	
				break;
			case 'p':
				ToolsFlags.flag_pattern = 1;
				ToolsFlags.pattern = strdup(optarg);
				break;
			case 'v':	
				ToolsFlags.flag_hide = 1;
				break;
			case 'n':	
				ToolsFlags.flag_columnName = 1;
				break;
			case 'l':	
				ToolsFlags.flag_simple = 1;
				break;
			case 'a':	
				ToolsFlags.flag_origin = 1;
				break;
			case 'd':	
				ToolsFlags.flag_dup = 1;
				break;
			case 'M':	
				ToolsFlags.flag_mapq = 1;
				break;
			case 'u':
				ToolsFlags.flag_coverage = 1;
				ToolsFlags.coverage_threshold = strdup(optarg);
				break;
			case 'w':
				ToolsFlags.flag_flow	= 1;
				break;
			case 'z':
				ToolsFlags.flag_zm	= 1;
				break;
			case 'h':	
				usage();
				return 0;
			case 'x':
				break;
			default: 	
				printf("Error: Parameter\n");	
				return -1;
		}
	}
	//Piroity 
	//Bed File > Region > All
	
	if (flag_bam == 1){	
		file_bam_i = fopen(bam,"rb");	
		if (access(bai_1, R_OK) != -1){	
			file_bai_i = fopen(bai_1,"rb");	
			//printf("bai_1:\t%s\n",bai_1);
		} else if (access(bai_2, R_OK) != -1){			
			file_bai_i = fopen(bai_2,"rb");	
			//printf("bai_2:\t%s\n",bai_2);
		}else {
			usage();
			return -1;
		}
	} else{	
		usage();
		return -1;	
	}
	if (flag_fasta == 1){
		file_fasta_i	= fopen(fasta, "r");	
		file_fai_i	= fopen(fai, "r");	
	}



	if (ToolsFlags.flag_bed == 1){	
		file_bed_i = fopen(bed,"r");
	}

	if (ToolsFlags.flag_target == 1){
		ToolsFlags.file_target = fopen(target, "r");
	}

	if (ToolsFlags.flag_start == 1){
		if (ToolsFlags.chromosome == NULL){
			printf("\n[Error!]\tChromosome is missed.\n");
			return -1;
		}
	}

	if (ToolsFlags.flag_end == 1){
		if (ToolsFlags.chromosome == NULL){
			printf("\n[Error!]\tChromosome is missed.\n");
			return -1;
		}
 		if (ToolsFlags.flag_start == 0){
			printf("\n[Error!]\tStart position is missed.\n");
			return -1;
		}
	}


	//start--;
	if (strcmp ( mode, "depthdist")==0){		BamDepthDist(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "indel")==0){
	}else if (strcmp ( mode, "ampsummary")==0){	BamAmp(file_bam_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "quality")==0){	BamSinglePointQuality(file_bam_i, file_bai_i, &ToolsFlags);
	}else if (strcmp ( mode, "stat")==0){		BamStat(file_bam_i, &ToolsFlags);
	}else if (strcmp ( mode, "pattern")==0){	BamPattern(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "ins")==0){		BamPoly(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "length")==0){		BamMappingLength(file_bam_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "del")==0){		BamDelTxt(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "purex")==0){		BamPureX(file_bam_i, file_bai_i, file_bed_i, &ToolsFlags);
	}else if (strcmp ( mode, "trim")==0){		
		if (ToolsFlags.flag_target == 0){	usage(); return -1;}
		BamTrim(file_bam_i, file_bed_i, &ToolsFlags);
//	}else if (strcmp ( mode, "iontorrent")==0){	
//		if (flag_fasta == 0){	usage(); return -1;}
//		if (ToolsFlags.flag_zm == 0 && ToolsFlags.flag_flow == 0){	usage(); return -1;}
//		ToolsFlags.flag_header	= 1;	BamTest(file_bam_i, file_bai_i, file_bed_i, file_fasta_i, file_fai_i, &ToolsFlags);
	}else {
		printf("[Warning]No mode or mode is error\n");
		usage();
		return -1;
	}


	if (flag_bam == 1){	fclose(file_bam_i);	fclose(file_bai_i);}
	if (flag_fasta == 1){	fclose(file_fasta_i);	fclose(file_fai_i);}
	if (ToolsFlags.flag_bed == 1){	fclose(file_bed_i);	}
//
	return 0;

}
