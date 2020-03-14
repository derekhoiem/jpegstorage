#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jpeglib.h>
#include <gmp.h>

const char IDENTIFIER[] = "XJ\0\1";
float BITS_PER_BLOCK_THRESHOLD;
const unsigned long int BLOCKS_PER_SECTION = 64;
const unsigned int SECTION_WARNING_THRESHOLD = 16*64;

float BASIS[64][64];
unsigned int STATES[64];
float STATES_LOG[64];
float REMAINING_BITS[64];
unsigned int ORDER[64];

int print_notice();
int read_quantization_table(const char file_name[], JQUANT_TBL * table);
int setup_bits(const char file_name[], unsigned long int * header_size, unsigned char ** header_bits, unsigned long int * raw_size, unsigned char ** raw_bits);
int write_jpeg(const char file_name[], const JQUANT_TBL * table, int ** blocks, const unsigned int number_of_blocks, const int image_width, const int image_height);
void initialize_tables(const JQUANT_TBL * table);
void initialize_basis(const JQUANT_TBL * table);
void dct_matrix(double matrix[64]);

int encode(mpz_t raw_bits, int ** blocks, const unsigned long int number_of_blocks);

void display_int_matrix(const int matrix[64]);
void display_float_matrix(const float matrix[64]);
void display_double_matrix(const double matrix[64]);
void display_jblock(const JBLOCK block);
void display_jquant_tbl(const JQUANT_TBL * table);

unsigned long int count_loops = 0;

/*
/  parameters:
/   example jpeg,
/   input file,
/   output jpeg
/   output aspect ratio
*/
const char example_file_name[] = "n1910415_24769.jpg";
const char input_file_name[] = "music.rar";
const char jpeg_file_name[] = "test.jpg";
const char aspect_ratio[] = "200x600";

int main()
{
	JQUANT_TBL quantization_table;
	unsigned char * raw_bits;
	unsigned long int raw_size;
	unsigned char * header_bits;
	unsigned long int header_size;

	unsigned int image_width = 200;
	unsigned int image_height = 600;

	unsigned long int number_of_blocks;
	unsigned long int number_of_sections;
	unsigned long int section_size;

	unsigned long int current_block;
	unsigned long int current_section;
	unsigned long int current_byte;

	int ** blocks;
	mpz_t gmp_bits;

	/* print notice */
	print_notice();

	/* read example quantization table */
	read_quantization_table(example_file_name, &quantization_table);

	/* initialize tables */
	initialize_tables(&quantization_table);

	/* read file and setup bits */
	setup_bits(input_file_name, &header_size, &header_bits, &raw_size, &raw_bits);

	/* calculate number of blocks, number of sections, and size of sections in image */
	number_of_blocks = image_width/8 * image_height/8;
	number_of_sections = (number_of_blocks + BLOCKS_PER_SECTION - 2) / BLOCKS_PER_SECTION;
	section_size = (raw_size * BLOCKS_PER_SECTION + number_of_blocks - 2) / (number_of_blocks - 1);
	BITS_PER_BLOCK_THRESHOLD = 8.0 * section_size / BLOCKS_PER_SECTION;

	/* do checks on block vs. bits basis */
	if (section_size > SECTION_WARNING_THRESHOLD)
		printf("Warning! The parameters given may be unachieveable.\n");

	/* allocate big integer */
	mpz_init2(gmp_bits, section_size);

	/* allocate all blocks */
	blocks = malloc(number_of_blocks*sizeof(int *));

	if (blocks == NULL)
		exit(1);

	for (current_block = 0; current_block < number_of_blocks; current_block++)
	{
		blocks[current_block] = malloc(64*sizeof(int));

		if (blocks == NULL)
			exit(1);

		memset(blocks[current_block], '\0', 64*sizeof(int));
	}

	/* encode raw bits in blocks */
	printf("Encoding %lu bytes in %lu blocks and %lu sections (%lu bytes per section)\n", header_size + raw_size, number_of_blocks, number_of_sections, section_size);

	/* encode header block */
	//printf("from: header, size: %lu, at: blocks[0], size: 1\n", header_size);

	mpz_import(gmp_bits, header_size, -1, 1, 0, 0, header_bits);

	if (encode(gmp_bits, &blocks[0], 1))
		printf("Encoded header successfully!\n");
	else
	{
		printf("Failed to encode header!\n");
		exit(1);
	}

	/* encode short section */
	//printf("from: raw_bits[0], size: %lu, at: blocks[1], size: %lu\n", raw_size % section_size, (number_of_blocks - 1) % BLOCKS_PER_SECTION);

	mpz_import(gmp_bits, raw_size % section_size, -1, 1, 0, 0, &raw_bits[0]);

	if (encode(gmp_bits, &blocks[1], (number_of_blocks - 1) % BLOCKS_PER_SECTION))
		printf("Encoded section 1 successfully!\n");
	else
	{
		printf("Failed to encode section 1!\n");
		exit(1);
	}

	current_block = (number_of_blocks - 1) % BLOCKS_PER_SECTION + 1;
	current_byte = raw_size % section_size;

	for (current_section = 0; current_section < number_of_sections - 1; current_section++)
	{
		//printf("from: raw_bits[%lu], size: %lu, at: blocks[%lu], size: %lu\n", current_byte, section_size, current_block, BLOCKS_PER_SECTION);

		mpz_import(gmp_bits, section_size, -1, 1, 0, 0, &raw_bits[current_byte]);
		//printf("%u %lu\n", mpz_sizeinbase(gmp_bits, 2), section_size*8);

		/* encode this section */
		if (encode(gmp_bits, &blocks[current_block], BLOCKS_PER_SECTION))
			printf("Encoded section %lu successfully!\n", current_section + 2);
		else
		{
			printf("Failed to encode section %lu!\n", current_section + 2);
			exit(1);
		}

		current_block += BLOCKS_PER_SECTION;
		current_byte += section_size;
	}

	/* free structures we are finished with */
	mpz_clear(gmp_bits);
	free(header_bits);
	free(raw_bits);

	//printf("Loops: %lu\n", count_loops);

	write_jpeg(jpeg_file_name, &quantization_table, blocks, number_of_blocks, image_width, image_height);

	/* free blocks */
	for (current_block = 0; current_block < number_of_blocks; current_block++)
		free(blocks[current_block]);

	free(blocks);

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

int read_quantization_table(const char file_name[], JQUANT_TBL * table)
{
	FILE * in_file;
	struct jpeg_decompress_struct jpeg_info;
	struct jpeg_error_mgr jpeg_error;

	/* open file */
	in_file = fopen(file_name, "rb");

	/* detect errors */
	if (in_file == NULL)
	{
		fprintf(stderr, "An error occured while opening \"%s\" for reading!\n", file_name);
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
	if (jpeg_info.quant_tbl_ptrs == NULL || jpeg_info.quant_tbl_ptrs[0] == NULL)
	{
		fprintf(stderr, "An error occured as the image quantization table is missing in \"%s\"\n", file_name);
		exit(1);
	}

	/* copy quantization table */
	memcpy(table, jpeg_info.quant_tbl_ptrs[0], sizeof(JQUANT_TBL));

	/* close file */
	fclose(in_file);

	/* free allocated memory */
	jpeg_destroy_decompress(&jpeg_info);

	return 0;
}

int setup_bits(const char file_name[], unsigned long int * header_size, unsigned char ** header_bits, unsigned long int * raw_size, unsigned char ** raw_bits)
{
	FILE * in_file;
	unsigned long int file_size;
	unsigned long int file_name_length;

	/* open file */
	in_file = fopen(file_name, "rb");

	/* detects errors */
	if (in_file == NULL)
	{
		printf("An error occured while opening \"%s\" for reading!\n", file_name);
		exit(1);
	}

	/* determine file size */
	fseek(in_file, 0, SEEK_END);
	file_size = ftell(in_file);
	fseek(in_file, 0, SEEK_SET);

	/* allocate array for raw bits */
	file_name_length = strlen(file_name);
	*header_size = 4 + sizeof(unsigned char) + sizeof(unsigned long int);
	*raw_size = file_name_length + file_size;
	*header_bits = malloc(*header_size);
	*raw_bits = malloc(*raw_size);

	if (*header_bits == NULL || *raw_bits == NULL)
	{
		printf("An error occured while attempting to allocate a memory buffer!\n");
		exit(1);
	}

	/* setup header bits */
	memcpy(*header_bits, IDENTIFIER, 4);
	(*header_bits)[4] = file_name_length;
	memcpy(*header_bits + 4 + sizeof(unsigned char), &file_size, sizeof(unsigned long int)); /* little endian */

	/* setup raw bits */
	memcpy(*raw_bits, file_name, file_name_length);
	fread(*raw_bits + file_name_length, sizeof(unsigned char), file_size, in_file);

	/* detect errors */
	if (!ferror(in_file) && !feof(in_file) && getc(in_file) == EOF)
	{
		printf("\"%s\" was read successfully!\n", file_name);
	}
	else
	{
		printf("An error occured while reading \"%s\"!\n", file_name);
		exit(1);
	}

	/* close file */
	fclose(in_file);

	return 0;
}

int write_jpeg(const char file_name[], const JQUANT_TBL * table, int ** blocks, const unsigned int number_of_blocks, const int image_width, const int image_height)
{
	FILE * out_file;
	struct jpeg_compress_struct jpeg_info;
	struct jpeg_error_mgr jpeg_error;

	jvirt_barray_ptr coefficient_blocks;
	JBLOCKARRAY row_blocks;

	unsigned long int current_block;
	unsigned int current_row;
	unsigned int current_column;
	unsigned int current_coefficient;


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
	coefficient_blocks = jpeg_info.mem->request_virt_barray((j_common_ptr)&jpeg_info, JPOOL_IMAGE, FALSE, image_width/8, image_height/8, 1);
	jpeg_info.mem->realize_virt_arrays((j_common_ptr)&jpeg_info);

	current_block = 0;

	for (current_row = 0; current_row < image_height/8; current_row++)
	{
		row_blocks = jpeg_info.mem->access_virt_barray((j_common_ptr)&jpeg_info, coefficient_blocks, current_row, 1, TRUE);

		for (current_column = 0; current_column < image_width/8; current_column++)
		{
			for (current_coefficient = 0; current_coefficient < 64; current_coefficient++)
			{
				row_blocks[0][current_column][ORDER[current_coefficient]] = blocks[current_block][current_coefficient];
			}

			current_block++;
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
	float log_sum;
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

	initialize_basis(table);

	log_sum = 0.0;
	for (i = 64 - 1; i >= 0; i--)
	{
		STATES_LOG[i] = log((float)STATES[i]) / log(2.0);
		log_sum += STATES_LOG[i];
		REMAINING_BITS[i] = log_sum;
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

int encode(mpz_t raw_bits, int ** blocks, const unsigned long int number_of_blocks)
{
	mpz_t branch_lefts[64];
	float branch_samples[64][64];
	float branch_bit_counts[64];

	mpz_t working_left;
	float working_samples[64];
	float working_bit_count;
	int working_block[64];
	int working_bits;

	float test_samples[64];

	int last_branches[64];
	int last_branch;

	mpz_t best_left;
	float best_bit_count;

	unsigned long int current_block;
	int current_coefficient;
	int i;
	int done;
	int range;

	/* initialize blocks to 0 */
	for (current_block = 0; current_block < number_of_blocks; current_block++)
		memset(blocks[current_block], '\0', 64*sizeof(int));

	/* start without success */
	done = 0;

	/* initialize large integers */
	mpz_init_set(best_left, raw_bits);
	mpz_init2(working_left, mpz_sizeinbase(raw_bits, 2));

	for (current_coefficient = 0; current_coefficient < 64; current_coefficient++)
		mpz_init2(branch_lefts[current_coefficient], mpz_sizeinbase(raw_bits, 2));

	for (current_block = 0; current_block < number_of_blocks && !done; current_block++)
	{
		/* starting new block */
		best_bit_count = 0.0;
		last_branch = -1;

		/* reset working entities */
		mpz_swap(best_left, working_left);
		memset(working_samples, '\0', 64*sizeof(float));
		working_bit_count = 0.0;
		memset(working_block, '\0', 64*sizeof(int));

		for (current_coefficient = 0; ; current_coefficient++)
		{
			count_loops++;

			/* check for running out of coefficients */
			if (current_coefficient >= 64 || (working_bit_count + REMAINING_BITS[current_coefficient]) < best_bit_count)
			{
				/*
				printf("Up:   ");
				mpz_out_str(stdout, 16, working_left);
				*/

				/* check for new best */
				if (working_bit_count > best_bit_count)
				{
					memcpy(blocks[current_block], working_block, 64*sizeof(int));
					best_bit_count = working_bit_count;
					mpz_swap(best_left, working_left);

					/* check if adequate */
					if (best_bit_count > BITS_PER_BLOCK_THRESHOLD)
					break;
				}

				/* check for finish */
				if (last_branch < 0)
					break;

				/* restore values from branch */
				mpz_swap(working_left, branch_lefts[last_branch]);
				memcpy(working_samples, branch_samples[last_branch], 64*sizeof(float));
				working_bit_count = branch_bit_counts[last_branch];
				working_block[last_branch] = 0;

				/* backtrack to one after last branch */
				current_coefficient = last_branch;
				last_branch = last_branches[last_branch];

				/*
				printf("->");
				mpz_out_str(stdout, 16, working_left);
				printf("\n");
				//display_int_matrix(working_block);
				*/
				continue;
			}

			/* may want to use: mp_limb_t mpn_divmod_1 (mp limb t *r1p, mp limb t *s2p, mp size t s2n, mp limb t s3limb) */
			/* determine bits this coefficient can encode */
			working_bits = (int)mpz_tdiv_q_ui(branch_lefts[current_coefficient], working_left, STATES[current_coefficient]);

			/* apply coefficient mapping */
			if (working_bits & 0x1)
				working_bits = -(working_bits + 1) / 2;
			else
				working_bits = (working_bits + 1) / 2 + 1;

			range = 1;

			/* add this number of the basis into samples so far while also checking range */
			for (i = 0; i < 64; i++)
			{
				test_samples[i] = working_samples[i] + BASIS[current_coefficient][i]*working_bits;

				/* check range, 0.5 rounds to nearest (this could be important when reversibility is verified) */
				if ((test_samples[i] + 0.5) < -128 || (test_samples[i] - 0.5) > 127)
				{
					range = 0;
					break;
				}
			}

			/* if not in range, move to next coefficient */
			if (!range)
				continue;

			/* will it finish */
			if (mpz_cmp_ui(branch_lefts[current_coefficient], 0) == 0)
			{
				/* finished encoding all bits, no need to worry about best */
				working_block[current_coefficient] = working_bits;
				memcpy(blocks[current_block], working_block, 64*sizeof(int));
				done = 1;
				break;
			}

			/* try this coefficient */
			mpz_swap(branch_lefts[current_coefficient], working_left);
			memcpy(branch_samples[current_coefficient], working_samples, 64*sizeof(float));
			branch_bit_counts[current_coefficient] = working_bit_count;
			last_branches[current_coefficient] = last_branch;

			memcpy(working_samples, test_samples, 64*sizeof(float));
			working_bit_count += STATES_LOG[current_coefficient];
			last_branch = current_coefficient;
			working_block[current_coefficient] = working_bits;

			/*
			printf("Down: ");
			mpz_out_str(stdout, 16, branch_lefts[current_coefficient]);
			printf("->");
			mpz_out_str(stdout, 16, working_left);
			printf(" (%u->%u)\n", last_branches[current_coefficient], last_branch);
			display_int_matrix(working_block);
			*/
		}
	}

	/* free anything allocated */
	for (current_coefficient = 0; current_coefficient < 64; current_coefficient++)
		mpz_clear(branch_lefts[current_coefficient]);

	mpz_clear(working_left);
	mpz_clear(best_left);

	return done;
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
