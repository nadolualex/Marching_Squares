
#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define Min(a, b) ((a) < (b) ? (a) : (b))

struct Thread {
    int thread_id;
    char *input_file;
    char *output_file;
    int threads_no;
    ppm_image **contour_map;
    ppm_image *scaled_image;
    unsigned char **grid;
    pthread_barrier_t *barrier;
};


// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }

    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }

    free(grid);

    free(image->data);
    free(image);
}

void *marching_squares (void *arg) {
    char **argv = (char **)arg;
    struct Thread *thread = (struct Thread *)arg;
    int thread_id = thread->thread_id;
    int threads_no = thread->threads_no;
    ppm_image **contour_map = thread->contour_map;
    ppm_image *scaled_image = thread->scaled_image;
    unsigned char **grid = thread->grid;
    pthread_barrier_t *barrier = thread->barrier;

    // Initialization
    int start = thread_id * CONTOUR_CONFIG_COUNT / threads_no;
    int end = Min((thread_id + 1) * CONTOUR_CONFIG_COUNT / threads_no, CONTOUR_CONFIG_COUNT);

    for (int i = start; i < end; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        contour_map[i] = read_ppm(filename);
    }

    ppm_image *image = read_ppm(argv[1]);

    pthread_barrier_wait(barrier);

    // Rescale image
    uint8_t sample[3];

    // We only rescale downwards
    if (!(image->x <= RESCALE_X && image->y <= RESCALE_Y)) {
        scaled_image->x = RESCALE_X;
        scaled_image->y = RESCALE_Y;

        scaled_image->data = (ppm_pixel*)malloc(scaled_image->x * scaled_image->y * sizeof(ppm_pixel));
        if (!scaled_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        // Use bicubic interpolation for scaling
        int start = thread_id * scaled_image->x / threads_no;
        int end = Min((thread_id + 1) * scaled_image->x / threads_no, scaled_image->x);
        for (int i = start; i < end; i++) {
            pthread_barrier_wait(barrier);
            for (int j = 0; j < scaled_image->y; j++) {
                float u = (float)i / (float)(scaled_image->x - 1);
                float v = (float)j / (float)(scaled_image->y - 1);
                sample_bicubic(image, u, v, sample);

                scaled_image->data[i * scaled_image->y + j].red = sample[0];
                scaled_image->data[i * scaled_image->y + j].green = sample[1];
                scaled_image->data[i * scaled_image->y + j].blue = sample[2];
            }
        }

        free(image->data);
        free(image);
    }
    else {
        scaled_image = image;
    }

    pthread_barrier_wait(barrier);

    // Grid
    int step_x = STEP;
    int step_y = STEP;

    int p = scaled_image->x / step_x;
    int q = scaled_image->y / step_y;
    
    start = thread_id * (p + 1) / threads_no;
    end = Min((thread_id + 1) * (p + 1) / threads_no, p + 1);

    for (int i = start; i < end; i++) {
        grid[i] = (unsigned char *)malloc((scaled_image->y / STEP + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    pthread_barrier_wait(barrier);

    start = thread_id * p / threads_no;
    end = Min((thread_id + 1) * p / threads_no, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = scaled_image->data[i * step_x * scaled_image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
        pthread_barrier_wait(barrier);
    }

    start = thread_id * q / threads_no;
    end = Min((thread_id + 1) * q / threads_no, q);

    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = scaled_image->data[i * step_x * scaled_image->y + scaled_image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }

    pthread_barrier_wait(barrier);

    start = thread_id * p / threads_no;
    end = Min((thread_id + 1) * p / threads_no, p);

    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = scaled_image->data[(scaled_image->x - 1) * scaled_image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    pthread_barrier_wait(barrier);


    // March
    if(scaled_image == image) {
        start = 0;
        end = p;
    } else {
        start = thread_id * p / threads_no;
        end = Min((thread_id + 1) * p / threads_no, p);
    }

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(scaled_image, contour_map[k], i * step_x, j * step_y);
        }
    }

    // Writing output file

    FILE *fp = fopen(argv[2], "wb");
    if (!fp) {
        fprintf(stderr, "Unable to open file '%s'\n", argv[2]);
        exit(1);
    }

    // Write the header file image format
    fprintf(fp, "P6\n");

    // image size
    fprintf(fp, "%d %d\n", scaled_image->x, scaled_image->y);

    // RGB component depth
    fprintf(fp, "%d\n", RGB_COMPONENT_COLOR);

    // pixel data
    fwrite(scaled_image->data, 3 * scaled_image->x, scaled_image->y, fp);
    fclose(fp);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    
    int threads_no = atoi(argv[3]);
    pthread_t tid[threads_no];

    if (argc < 3) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file>\n");
        return 1;
    }

    ppm_image **contour_map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!contour_map) {
        fprintf(stderr, "Unable to allocate memory\n");
        return 1;
    }

    ppm_image *scaled_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!scaled_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        return 1;
    }

    unsigned char **grid = (unsigned char **)malloc((scaled_image->x / STEP + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        return 1;
    }

    struct Thread thread[threads_no];
    
    pthread_barrier_t barrier;

    pthread_barrier_init(&barrier, NULL, threads_no);

    for(int i = 0; i < threads_no; i ++) {
        thread[i].input_file = argv[1];
        thread[i].output_file = argv[2];
        thread[i].thread_id = i;
        thread[i].threads_no = threads_no;
        thread[i].contour_map = contour_map;
        thread[i].scaled_image = scaled_image;
        thread[i].grid = grid;
        thread[i].barrier = &barrier;

        pthread_create(&tid[i], NULL, marching_squares, &thread[i]);
    }

    for(int i = 0; i < threads_no; i ++) {
        pthread_join(tid[i], NULL);
    }

    pthread_barrier_destroy(&barrier);

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for(int i = 0; i < scaled_image->x / STEP + 1; i ++) {
        free(grid[i]);
    }

    return 0;
}


