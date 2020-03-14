#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jpeglib.h>
#include <gmp.h>

const unsigned long int BLOCKS_PER_SECTION = 64;
const char IDENTIFIER[] = "XJ\0\1";

float BASIS[64][64];
unsigned int STATES[64];
unsigned int ORDER[64];

int print_notice();
void initialize_tables(const JQUANT_TBL * table);
void initialize_basis(const JQUANT_TBL * table);
void dct_matrix(double matrix[64]);

int decode(mpz_t raw_bits, const JBLOCKROW blocks, const unsigned int number_of_blocks);

void display_int_matrix(const int matrix[64]);
void display_float_matrix(const float matrix[64]);
void display_double_matrix(const double matrix[64]);
void display_jblock(const JBLOCK block);
void display_jquant_tbl(const JQUANT_TBL * table);

int main(int argc, char * argv[])
{
	/*
	/  parameters:
	/   input file
	*/

	/* input file */
	FILE * in_file;
	FILE * out_file;
	char * file_name;
	char * out_file_name;
	unsigned char file_name_length;
	unsigned int file_size;

	struct jpeg_decompress_struct jpeg_info;
	struct jpeg_error_mgr jpeg_error;
	jvirt_barray_ptr * jpeg_coef_arrays;
	JBLOCKARRAY blocks;
	unsigned char header[9];

	unsigned char * raw_bits;
	unsigned long int raw_size;
	mpz_t gmp_bits;

	unsigned long int number_of_blocks;
	long int blocks_left;
	unsigned int current_row;
	unsigned int current_column;
	unsigned long int blocks_on_row;

	unsigned long int number_of_sections;
	unsigned long int current_section;
	unsigned long int section_size;

	unsigned int width_in_blocks;
	unsigned int height_in_blocks;

	unsigned long int current_byte;

	/* print notice */
	print_notice();

	/* parse parameters */
	if (argc != 2
	 || strcmp(argv[1], "--help") == 0
	 || strcmp(argv[1], "--h") == 0
	 || strcmp(argv[1], "-help") == 0
	 || strcmp(argv[1], "-h") == 0
	 || strcmp(argv[1], "/?") == 0
	 || strcmp(argv[1], "/h") == 0
	 || strcmp(argv[1], "/help") == 0)
	{
		printf("Usage: %s <file name>.jpg\n", argv[0]);
		exit(0);
	}

	/* redirect file name */
	file_name = argv[1];

	/* open file */
	in_file = fopen(file_name, "rb");

	/* detect errors */
	if (in_file == NULL)
	{
		fprintf(stderr, "Error! Could not open \"%s\" for reading!\n", file_name);
		exit(1);
	}

	/* allocate and initialize JPEG decompression object */
	jpeg_info.err = jpeg_std_error(&jpeg_error);
	jpeg_create_decompress(&jpeg_info);

	/* specify data source */
	jpeg_stdio_src(&jpeg_info, in_file);

	/* read file parameters */
	jpeg_read_header(&jpeg_info, TRUE);

	/* detect errors */
	if (jpeg_info.image_width % 8 != 0 || jpeg_info.image_height % 8 != 0)
	{
		fprintf(stderr, "Error! Image dimensions are currently unsupported!\n");
		exit(1);
	}
	if (jpeg_info.quant_tbl_ptrs == NULL || jpeg_info.quant_tbl_ptrs[0] == NULL)
	{
		fprintf(stderr, "Error! Image quantization table is missing!\n");
		exit(1);
	}

	/* initialize tables */
	initialize_tables(jpeg_info.quant_tbl_ptrs[0]);

	/* read raw coefficients */
	jpeg_coef_arrays = jpeg_read_coefficients(&jpeg_info);

	/* read header */
	blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, jpeg_coef_arrays[0], 0, 1, FALSE);

	//display_jblock(blocks[0][0]);

	/* initialize bits for header */
	memset(header, '\0', 9);
	mpz_init2(gmp_bits, 9);

	/* decode header */
	mpz_set_ui(gmp_bits, 0);
	decode(gmp_bits, &blocks[0][0], 1);

	/* detect errors */
	if (mpz_sizeinbase(gmp_bits, 2) > 9*8)
	{
		fprintf(stderr, "Error! Could not read header section!\n");
		exit(1);
	}

	/* copy bits */
	mpz_export(header, NULL, -1, 1, 0, 0, gmp_bits);

	/* ensure correct identifier */
	if (memcmp(header, IDENTIFIER, 4) != 0)
	{
		fprintf(stderr, "Error! Could not identify the encoding format\n");
		exit(1);
	}

	/* load sizes */
	file_name_length = header[4];
	memcpy(&file_size, &header[5], sizeof(unsigned int));

	//printf("\n%d\n%lu", file_name_length, file_size);

	/* allocate array for raw bits */
	raw_size = file_name_length + file_size;
	raw_bits = malloc(raw_size);
	out_file_name = malloc(file_name_length + 1);

	if (raw_bits == NULL || out_file_name == NULL)
	{
		fprintf(stderr, "Error! Could not allocate memory buffers!\n");
		exit(1);
	}

	/* initialize raw bits */
	out_file_name[file_name_length] = '\0';
	memset(raw_bits, '\0', raw_size);

	/* calculate number of blocks, number of sections, and size of sections in image */
	width_in_blocks = jpeg_info.comp_info[0].width_in_blocks;
	height_in_blocks = jpeg_info.comp_info[0].height_in_blocks;
	number_of_blocks = width_in_blocks * height_in_blocks;
	number_of_sections = (number_of_blocks + BLOCKS_PER_SECTION - 2) / BLOCKS_PER_SECTION;
	section_size = (raw_size * BLOCKS_PER_SECTION + number_of_blocks - 2) / (number_of_blocks - 1);

	printf("Decoding %lu bytes in %lu blocks and %lu sections (%lu bytes per section)\n", raw_size, number_of_blocks, number_of_sections, section_size);

	/* reallocate big integer */
	mpz_realloc2(gmp_bits, section_size);

	current_row = height_in_blocks - 1;
	current_column = 0;
	current_byte = file_size - section_size + file_name_length;
	blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, jpeg_coef_arrays[0], current_row, 1, FALSE);

	/* decode full sections */
	for (current_section = number_of_sections - 1; current_section > 0; current_section--)
	{
		mpz_set_ui(gmp_bits, 0);
		blocks_left = BLOCKS_PER_SECTION;

		do
		{
			//printf(" r: %u", current_row);

			blocks_on_row = width_in_blocks - current_column;

			if (blocks_on_row > blocks_left)
				blocks_on_row = blocks_left;

			//printf(" b: %lu", blocks_on_row);
			//printf(" c: %u", current_column);
			//printf(" c: %u", width_in_blocks - current_column - 1);
			//printf(" d: %lu\n", width_in_blocks - current_column - blocks_on_row);

			decode(gmp_bits, &blocks[0][width_in_blocks - current_column - blocks_on_row], blocks_on_row);

			if ((current_column + blocks_on_row) == width_in_blocks)
			{
				current_row--;
				blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, jpeg_coef_arrays[0], current_row, 1, FALSE);
			}

			current_column += blocks_on_row;
			current_column %= width_in_blocks;

			blocks_left -= blocks_on_row;
		} while (blocks_left > 0);

		//printf("%u %lu\n", mpz_sizeinbase(gmp_bits, 2), section_size*8);

		/* detect errors */
		if (mpz_sizeinbase(gmp_bits, 2) - 8 > section_size*8)
		{
			fprintf(stderr, "Error! Could not read section %lu!\n", current_section + 1);
			exit(1);
		}

		//printf("current byte: %lu\n", current_byte);
		/* write section */
		mpz_export(&raw_bits[current_byte], NULL, -1, 1, 0, 0, gmp_bits);
		current_byte -= section_size;
	}

	/* decode short section */
	mpz_set_ui(gmp_bits, 0);
	blocks_left = (number_of_blocks - 1) % BLOCKS_PER_SECTION;

	do
	{
		blocks_on_row = width_in_blocks - current_column;

		if (blocks_on_row > blocks_left)
			blocks_on_row = blocks_left;

		decode(gmp_bits, &blocks[0][width_in_blocks - current_column - blocks_on_row], blocks_on_row);

		if ((current_column + blocks_on_row) == width_in_blocks)
		{
			current_row--;
			blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, jpeg_coef_arrays[0], current_row, 1, FALSE);
		}

		current_column += blocks_on_row;
		current_column %= width_in_blocks;

		blocks_left -= blocks_on_row;
	} while (blocks_left > 0);

	/* detect errors */
	if (mpz_sizeinbase(gmp_bits, 2) > (raw_size % section_size)*8)
	{
		fprintf(stderr, "Error! Could not read section 1!\n");
		exit(1);
	}

	/* write section */
	mpz_export(&raw_bits[0], NULL, -1, 1, 0, 0, gmp_bits);

	/* get file name */
	memcpy(out_file_name, raw_bits, file_name_length);

	/* test if file already exists */
	out_file = fopen(out_file_name, "rb");

	if (out_file != NULL)
	{
		printf("Error! \"%s\" already exists!\n", out_file_name);
		fclose(out_file);
		exit(1);
	}

	/* open file */
	out_file = fopen(out_file_name, "wb");

	/* detect errors */
	if (out_file == NULL)
	{
		printf("Error! Could not open \"%s\" for writing!\n", out_file_name);
		exit(1);
	}

	printf("Writing \"%s\"\n", out_file_name);

	/* output file */
	fwrite(&raw_bits[file_name_length], raw_size - file_name_length, 1, out_file);

	/* close file */
	fclose(out_file);

	/* free structures we are finished with */
	mpz_clear(gmp_bits);
	free(raw_bits);
	free(out_file_name);
	jpeg_destroy_decompress(&jpeg_info);

	return 0;
}

int print_notice()
{
	printf("This software is based in part on the work of the Independent JPEG Group\n");
	printf("This software uses the GNU MP Bignum Library\n");
	printf("  as covered by the Lesser General Public License (LGPL).\n");
	printf("  You should have received a copy of the LGPL with this software.\n");
	printf("  It can also be found at http://www.gnu.org/copyleft/lesser.html\n");
	printf("  The LGPL grants you certain rights, including, but not limited to,\n");
	printf("  source code for the library and the unlinked object code of this program\n");
	printf("  The library can also be found at http://www.swox.com/gmp/\n");
	printf("\n");
	printf("\n");

	return 0;
}

int decode(mpz_t raw_bits, const JBLOCKROW blocks, const unsigned int number_of_blocks)
{
	long int current_block;
	int i;
	JCOEF current_coefficient;

	current_block = (signed)number_of_blocks;

	for (current_block--; current_block >= 0; current_block--)
	{
		for (i = 64 - 1; i >= 0; i--)
		{
			current_coefficient = blocks[current_block][ORDER[i]];
			if (current_coefficient != 0)
			{
				mpz_mul_ui(raw_bits, raw_bits, STATES[i]);
				if (current_coefficient < 0)
					mpz_add_ui(raw_bits, raw_bits, -current_coefficient*2 - 1);
				else
					mpz_add_ui(raw_bits, raw_bits, (current_coefficient - 1)*2);
			}
		}
	}

	return 1;
}

struct pair
{
	unsigned int states;
	unsigned int order;
};

int compare_pairs(const void * a, const void * b)
{
	return (((struct pair *)a)->states - ((struct pair *)b)->states);
}

void initialize_tables(const JQUANT_TBL * table)
{
	int i;
	int j;
	int max;
	struct pair pairs[64];

	for (i = 0; i < 64; i++)
	{
		STATES[i] = 2048 / table->quantval[i] - 1;
		ORDER[i] = i;
	}

	initialize_basis(table);

	for (i = 0; i < 64; i++)
	{
		max = 0;
		for (j = 1; j < 64; j++)
			if (BASIS[i][j] > BASIS[i][max])
				max = j;

		if (STATES[i] > 256 / BASIS[i][max])
			STATES[i] = 256 / BASIS[i][max];
	}

	for (i = 0; i < 64; i++)
	{
		pairs[i].states = STATES[i];
		pairs[i].order = i;
	}

	qsort(pairs, 64, sizeof(struct pair), compare_pairs);

	for (i = 0; i < 64; i++)
	{
		STATES[i] = pairs[i].states;
		ORDER[i] = pairs[i].order;
	}
}

void initialize_basis(const JQUANT_TBL * table)
{
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int m;
	unsigned int n;
	double dct[64];
	double intermediate;

	dct_matrix(dct);

	for (i = 0; i < 64; i++)
	{
		j = ORDER[i];
		k = 0;
		for (m = 0; m < 8; m++)
		{
			intermediate = dct[(j/8)*8 + m]*table->quantval[j];

			for (n = 0; n < 8; n++)
			{
				BASIS[i][k] = dct[j%8*8+n]*intermediate;
				k++;
			}
		}
	}
}

void dct_matrix(double matrix[64])
{
	unsigned int i;
	unsigned int j;
	unsigned int m;
	unsigned int n;
	double k;

	i = 0;
	for (m = 0; m < 8; m++)
	{
		if (m == 0)
			k = 0.35355339059327380;
		else
			k = 0.5;

		j = 1;
		for (n = 0; n < 8; n++)
		{
			matrix[i] = k*cos(0.19634954084936207*j*m);
			i++;
			j += 2;
		}
	}
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
