#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




//int	BamAmp(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags){
int	BamAmp(FILE *file_bam_i, FILE *file_bed_i, toolsFlags *ToolsFlags){

	FILE	*file_region_info_o;
	FILE	*file_region_stat_o;
//	FILE	*file_base_stat_o;
//	FILE	*file_base_info_o;

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t *address;
	uint8_t	*top;
	
	uint32_t	i;
	int	j;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;


	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;

	uint64_t	valid_region = 0;
	uint64_t	valid_all_region = 0;
	int	chr_in_target = 0;
	uint64_t	map_in = 0;
	uint64_t	map_out = 0;
	uint64_t	map_all_in = 0;
	uint64_t	map_all_out = 0;

	bamHeader	BamHeader;
	bedTable	BedTable;
	regionInformationLite	*RegionInformation;


	covBaseDistribution BaseDist_All;
	covRegionDistribution RegionDist_0_All;

	uint32_t	num_target = 0;	

//	int allocate_s_ptr	= 0;
	uint32_t allocate_e_ptr	= 0;	// the next one need to be allocate memory




/////////////////////////////////////////////////////
	//Start 
	////Bam Header
	////
	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536*2,sizeof(uint8_t));

	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/
	
	top = CatchBamHeader_v1 (file_bam_i, &BamHeader, stream_i, stream_o);
	counter	= top - stream_o;
	address = stream_o;

	/* comment on 2020/03/20 YC
	if (ToolsFlags->flag_hide == 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/

	//Bed File
	ToolsFlags->flag_unmerged	= 1;
	num_target = createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);
//	RegionInformation	= calloc(num_target, sizeof(regionInformation));
	RegionInformation	= calloc(num_target, sizeof(regionInformationLite));
	//printf("%d\n", num_target);
	//getchar();
/*
	for (i = 0; i < num_target; i++){
		RegionInformation[i].length	= BedTable.end[i] - BedTable.start[i];
		RegionInformation[i].ALL	= calloc(BedTable.end[i] - BedTable.start[i], sizeof(uint32_t));
		RegionInformation[i].DEL	= calloc(BedTable.end[i] - BedTable.start[i], sizeof(uint32_t));
		RegionInformation[i].overlap	= calloc(BedTable.end[i] - BedTable.start[i], sizeof(uint8_t));
		RegionInformation[i].start	= 0;
	}
*/	
//////////////////////////////////////////////////////////////////////////////////////////
	

	memset(&BaseDist_All, 0, sizeof(covBaseDistribution));
	memset(&RegionDist_0_All, 0, sizeof(covRegionDistribution));

	file_region_info_o 	= fopen("Region_info.txt","w");
	file_region_stat_o 	= fopen("Region_stat.txt","w");
//	file_base_stat_o 	= fopen("Base_stat.txt","w");
//	file_base_info_o 	= fopen("Base_info.txt","w");
	


	//Start Unzip and Process Bam File
	while ( (len_data = BGFZBlock(file_bam_i)) > 0){
		memmove( stream_o, address, sizeof(uint8_t)*counter);
		decompressBlock(&infstream, stream_i, stream_o+counter, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		
		top	= stream_o + counter + BgfzTail.I_size;
		address = stream_o;
		
		if (top - address > (intptr_t) sizeof(alignmentHeader)){
			memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
		}

		uint32_t	i_start = 0;
		while ( top - address >= (AlignmentHeader.block_size + 4) ){

			//Chromosome Change
			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (ToolsFlags->flag_hide == 1 ){
						printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
					}
					if (chr_in_target == 1){
/*						PrintBaseDist (file_base_stat_o, RegionInformation, &BedTable, BamHeader.chr_name[ref_ID], ref_ID, &BaseDist_All , 1);
						PrintRegionDist (file_stat_0_o  , file_ratio_0_o   ,&BedTable, RegionInformation, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 1);
*/						chr_in_target = 0;
						for (i = BedTable.table_start[ref_ID]; i < BedTable.table_end[ref_ID]; i++){
//							if (RegionInformation[i].count > 0){
//								free(RegionInformation[i].ALL);
//							}
						}
					}
					fprintf(file_region_stat_o,"%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);

					//printf("%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out, valid_region);
					map_all_in += map_in;
					map_all_out += map_out;
					valid_all_region += valid_region;
					valid_region	= 0;
					map_in	= 0;
					map_out	= 0;
				}	
				uint32_t max_end	= 0;

				if ( AlignmentHeader.refID != -1){	
					if (BedTable.table_start[AlignmentHeader.refID] != BedTable.table_end[AlignmentHeader.refID]){

						chr_in_target = 1;
						i_start	= BedTable.table_start[AlignmentHeader.refID];	
						// forward process
						for (i = BedTable.table_start[AlignmentHeader.refID]; i < BedTable.table_end[AlignmentHeader.refID];i++){
							if (BedTable.start[i] > max_end){
								valid_region	+= (BedTable.end[i] - BedTable.start[i]);
								max_end	= BedTable.end[i];
							}else if (BedTable.end[i] > max_end){
								valid_region	+= (BedTable.end[i] - max_end);
								max_end = BedTable.end[i];
							}
						}
						//forward process (second phase)
					}
				}
				ref_ID = AlignmentHeader.refID;	
				
			}
			if (ref_ID == -1){
				break;	
			}
			address += sizeof(alignmentHeader);	
			if ((AlignmentHeader.FLAG&4) == 0 && ((AlignmentHeader.FLAG&1024) == 0 || ToolsFlags->flag_dup==0)){
//			if ((AlignmentHeader.FLAG&4) == 0 && (AlignmentHeader.FLAG&256) == 0 && ((AlignmentHeader.FLAG&1024) == 0 || ToolsFlags->flag_dup==1)){
				
				if (chr_in_target == 1){
					uint32_t read_end	= alignment_EndPosition(address, &AlignmentHeader, BamHeader.chr_length[ref_ID]);
					int flag_ontarget	= 0;

					for (i = i_start; i < BedTable.table_end[AlignmentHeader.refID];i++){
						if ((read_end > BedTable.start[i]) && ((uint32_t) AlignmentHeader.pos < BedTable.end[i])){
							flag_ontarget	= 1;
							map_in++;
							break;
						}else if (BedTable.start[i] >= read_end ){
							break;
						}else if (BedTable.end[i] <= (uint32_t) AlignmentHeader.pos){
							CalculateStatInfo(file_region_info_o, &BedTable, RegionInformation, i_start, &BamHeader, ref_ID);
							if (RegionInformation[i_start].count > 0){
								free(RegionInformation[i_start].ALL);
							}
							i_start	= i+1;
						}
					}
					if (flag_ontarget == 0){
						map_out++;
					}
					if (flag_ontarget == 1){
						for (i = i_start; i < BedTable.table_end[AlignmentHeader.refID];i++){
							if (BedTable.start[i] >= read_end ){
								break;
							}
							
							if (i >= allocate_e_ptr){
								RegionInformation[i].ALL	= calloc(BedTable.end[i] - BedTable.start[i], sizeof(uint32_t));
								allocate_e_ptr	= i + 1;
							}
							int virtual_start	= (BedTable.start[i] < (uint32_t) AlignmentHeader.pos)? AlignmentHeader.pos - BedTable.start[i]: 0;
							int virtual_end	= (BedTable.end[i] < read_end) ? BedTable.end[i] - BedTable.start[i]: read_end - BedTable.start[i];
							for (j = virtual_start; j < virtual_end; j++){
								RegionInformation[i].ALL[j]++;
							}
							RegionInformation[i].count++;
						}
					}
				}else {
					map_out++;	
				}
			}
			address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

			if (top - address > (intptr_t) sizeof(alignmentHeader)){
				memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
			}else {
				break;	
			}
		}
		counter = top - address;
		if (ref_ID == -1){
			break;	
		}
	}
	
	if (ref_ID != -1){
		if (ToolsFlags->flag_hide == 1){	printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1,BamHeader.n_ref,BamHeader.chr_name[ref_ID]);	}

		if (chr_in_target == 1){
/*			PrintBaseDist (file_base_stat_o, RegionInformation, &BedTable, BamHeader.chr_name[ref_ID], ref_ID, &BaseDist_All , 1);
			PrintRegionDist (file_stat_0_o  , file_ratio_0_o   ,&BedTable, RegionInformation, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 1);
*/			chr_in_target = 0;
		}
		fprintf(file_region_stat_o,"%s\t%lu\t%lu\t%lu\n",BamHeader.chr_name[ref_ID],map_in,map_out,valid_region);
		map_all_in += map_in;
		map_all_out += map_out;
		valid_all_region += valid_region;
	}
		
/*	PrintBaseDist (file_base_stat_o, RegionInformation, &BedTable, "", 0, &BaseDist_All , 2);
	PrintRegionDist (file_stat_0_o  , file_ratio_0_o   ,&BedTable, RegionInformation, &BamHeader, ref_ID, 0  , &RegionDist_0_All  , 2);
*/	
	fprintf(file_region_stat_o,"Total\t%lu\t%lu\t%lu\t%f\t%f\n",
		map_all_in,
		map_all_out,
		valid_all_region,
		(float)map_all_in /(map_all_in+map_all_out),
		(float)map_all_out /(map_all_in+map_all_out));

	free(stream_i);	
	free(stream_o);	
	free(RegionInformation);	
	return 0;
}

