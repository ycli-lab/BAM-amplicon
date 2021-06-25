#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"

int	BamPureX(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags){


	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*buffer;

	int	j;

	int	len_data;
	int	ref_ID	= -1;


	z_stream infstream;
	bgfzTail BgfzTail;
	
	bamHeader	BamHeader;
	bedTable	BedTable;
	baiTable	BaiTable;
	uint64_t	offset_beg;
	uint64_t	offset_bgzf;


	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536,sizeof(uint8_t));
	buffer	= calloc(65536*2,sizeof(uint8_t));	

	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/
	
	if (ToolsFlags->flag_start == 1 || ToolsFlags->flag_end == 1 || ToolsFlags->flag_chromosome == 1 || ToolsFlags->flag_bed == 1){
		printf("This mode does not temporarily support the following parameters {'c', 's', 'e', 'r'}\n");
		return 0;
	}

	CatchBamHeader_v1 (file_bam_i, &BamHeader, stream_i, stream_o);

	BuildBaiTable (file_bai_i, &BaiTable);
	
	/* comment 2020/03/30 YC
	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/
	//Start Unzip and Process Bam File
	createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);
	for (j = 0; j < BamHeader.n_ref; j++){
		ref_ID	= j;
		if (BedTable.table_end[ref_ID] != BedTable.table_start[ref_ID]){
			offset_beg = ReturnOffset ( &BaiTable, ref_ID, BedTable.start[BedTable.table_start[ref_ID]]);
			offset_bgzf = offset_beg >> 16;

			if (offset_beg != 0){
				//Offset Bam File
				fseek(file_bam_i,offset_bgzf,SEEK_SET);
				len_data = BGFZBlock(file_bam_i);
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
//				fread(stream_i+2,sizeof(uint8_t),len_data,file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			}

			//Start Unzip and Process Bam File	
			while ( (len_data = BGFZBlock(file_bam_i)) > 0){
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
//				fread(stream_i+2,sizeof(uint8_t),len_data,file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail), 1, file_bam_i);
			}
		}
		break;
	}
	
	free(stream_i);	
	free(stream_o);	
	free(buffer);	
	return 0;
}


