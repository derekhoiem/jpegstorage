#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jpeglib.h>
#include <gmp.h>

int write_jpeg(const char file_name[], const JQUANT_TBL * table);

void display_int_matrix(const int matrix[64]);
void display_float_matrix(const float matrix[64]);
void display_double_matrix(const double matrix[64]);
void display_jblock(const JBLOCK block);
void display_jquant_tbl(const JQUANT_TBL * table);

const char jpeg_file_name[] = "test.jpg";

int main()
{
	JQUANT_TBL quantization_table;
	unsigned int current_coefficient;

	for (current_coefficient = 0; current_coefficient < 64; current_coefficient++)
		quantization_table.quantval[current_coefficient] = 8;


	write_jpeg(jpeg_file_name, &quantization_table);

	/* free blocks */

	return 0;
}

int write_jpeg(const char file_name[], const JQUANT_TBL * table)
{
	FILE * out_file;
	struct jpeg_compress_struct jpeg_info;
	struct jpeg_error_mgr jpeg_error;

	jvirt_barray_ptr coefficient_blocks;
	JBLOCKARRAY row_blocks;

	unsigned int current_row;
	unsigned int current_column;
	unsigned int image_width = 120;
	unsigned int image_height = 120;


	/* open file */
	out_file = fopen(file_name, "wb");

	/* detect errors */
	if (out_file == NULL)
	{
		fprintf(stderr, "An error occured while opening \"%s\" for writing!\n", file_name);
		exit(1);
	}

	/* allocate and initialize JPEG decompression object */
	jpeg_info.err = jpeg_std_error(&jpeg_error);
	jpeg_create_compress(&jpeg_info);

	/* specify data source */
	jpeg_stdio_dest(&jpeg_info, out_file);

	/* specify parameters */
	jpeg_info.image_width = image_width;
	jpeg_info.image_height = image_height;
	jpeg_info.input_components = 1;
	jpeg_info.in_color_space = JCS_GRAYSCALE;
	jpeg_info.num_components = 1;
	jpeg_set_defaults(&jpeg_info);

	if (jpeg_info.quant_tbl_ptrs == NULL || jpeg_info.quant_tbl_ptrs[0] == NULL)
	{
		fprintf(stderr, "An error occured while trying to select the quantization table!\n");
		exit(1);
	}

	memcpy(jpeg_info.quant_tbl_ptrs[0], table, sizeof(JQUANT_TBL));

	/* allocate a full-image virtual array */
	coefficient_blocks = jpeg_info.mem->request_virt_barray((j_common_ptr)&jpeg_info, JPOOL_IMAGE, TRUE, image_width/8, image_height/8, 1);
	jpeg_info.mem->realize_virt_arrays((j_common_ptr)&jpeg_info);

	for (current_row = 0; current_row < 8; current_row++)
	{
		row_blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, coefficient_blocks, current_row*2, 1, TRUE);

		for (current_column = 0; current_column < 8; current_column++)
		{
			row_blocks[0][current_column*2][current_row*8 + current_column] = 64;
		}

		if (current_row != 7)
		{
			row_blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, coefficient_blocks, current_row*2 + 1, 1, TRUE);
		}
	}

	/* output jpeg */
	jpeg_write_coefficients(&jpeg_info, &coefficient_blocks);
	jpeg_finish_compress(&jpeg_info);

	/* close file */
	fclose(out_file);

	/* free allocated memory */
	jpeg_destroy_compress(&jpeg_info);

	return 0;
}

void display_int_matrix(const int matrix[64])
{
	unsigned int i;

	for (i = 0; i < 64; i++)
	{
		printf("%9d", matrix[i]);

		if ((i + 1) % 8 == 0)
			printf("\n");
	}
	printf("\n\n");
}

void display_float_matrix(const float matrix[64])
{
	unsigned int i;

	for (i = 0; i < 64; i++)
	{
		printf("%9.5f", matrix[i]);

		if ((i + 1) % 8 == 0)
			printf("\n");
	}
	printf("\n\n");
}

void display_double_matrix(const double matrix[64])
{
	unsigned int i;

	for (i = 0; i < 64; i++)
	{
		printf("%9.5f", matrix[i]);

		if ((i + 1) % 8 == 0)
			printf("\n");
	}
	printf("\n\n");
}

void display_jblock(const JBLOCK block)
{
	unsigned int i;

	for (i = 0; i < 64; i++)
	{
		printf("%9d", block[i]);

		if ((i + 1) % 8 == 0)
			printf("\n");
	}
	printf("\n\n");
}

void display_jquant_tbl(const JQUANT_TBL * table)
{
	unsigned int i;

	for (i = 0; i < 64; i++)
	{
		printf("%9d", table->quantval[i]);

		if ((i + 1) % 8 == 0)
			printf("\n");
	}
	printf("\n\n");
}
