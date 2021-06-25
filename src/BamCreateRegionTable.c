#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BamCommonLibrary.h"

int	createRegionTable(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader, toolsFlags *ToolsFlags){

	int	num_target	= 0;
	int	j;

	int	start	= ToolsFlags->start;
	int	end	= ToolsFlags->end;
	char	*chromosome	= ToolsFlags->chromosome;

	if (file_bed_i != NULL){
		if (ToolsFlags->flag_unmerged){
			num_target = decodeBedFile(file_bed_i, BedTable, BamHeader);
		}else {
			num_target = decodeBedFile_Merge(file_bed_i, BedTable, BamHeader);
		}
	}else if (end > 0){
		//printf("Region\n");
		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		
		BedTable->start = calloc (1, sizeof(uint32_t));
		BedTable->end = calloc (1, sizeof(uint32_t));
		num_target = 1;
		BedTable->start[0] = start-1;
		BedTable->end[0] = end;
		
		for(j = 0; j < BamHeader->n_ref; j++){
			if (strcmp(chromosome, BamHeader->chr_name[j])==0){
				BedTable->table_start[j] = 0;
				BedTable->table_end[j] = 1;
				BedTable->table_max[j] = end;
				if( start > BamHeader->chr_length[j] || end > BamHeader->chr_length[j] || start > end){
					printf("[Error] Position is bigger than Size of %s\n",chromosome);
					exit(-1);	
				}
				break;
			}
		}
	}else if (start > 0){
		//printf("SinglePoint\n");
		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->start = calloc (1, sizeof(uint32_t));
		BedTable->end = calloc (1, sizeof(uint32_t));
		num_target = 1;
		BedTable->start[0] = start-1;
		BedTable->end[0] = start;
		
		for(j = 0; j < BamHeader->n_ref; j++){
			if (strcmp(chromosome, BamHeader->chr_name[j])==0){
				BedTable->table_start[j] = 0;
				BedTable->table_end[j] = 1;
				BedTable->table_max[j] = start+1;
				if( start > BamHeader->chr_length[j]){
					printf("[Error] Position is bigger than Size of %s\n",chromosome);
					exit(-1);	
				}
				break;
			}
		}	
	}else {
		//printf("Total\n");
		BedTable->table_start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->table_max = calloc (BamHeader->n_ref, sizeof(uint32_t));
		
		BedTable->start = calloc (BamHeader->n_ref, sizeof(uint32_t));
		BedTable->end = calloc (BamHeader->n_ref, sizeof(uint32_t));
		num_target = BamHeader->n_ref;
			
		for(j = 0; j < BamHeader->n_ref; j++){
			BedTable->start[j] = 0;
			BedTable->end[j] = BamHeader->chr_length[j];
			BedTable->table_start[j] = j;
			BedTable->table_end[j] = j+1;
			BedTable->table_max[j] = BamHeader->chr_length[j];
		}		
	}
	return num_target;
}
