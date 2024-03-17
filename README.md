# Parallel Marching Squares Algorithm with Pthreads

This repository contains an implementation of the Marching Squares algorithm, a technique used for extracting contours from images, optimized for parallel execution using Pthreads.

## Overview

The Marching Squares algorithm is a popular method for contour extraction in digital image processing. It identifies contours in a binary image by analyzing the values of neighboring pixels. While the algorithm is efficient, its computational demands can be significant, especially for large images.

To address this challenge, this project introduces a parallelized version of the Marching Squares algorithm using Pthreads, enabling faster processing of images through concurrent execution.

## Features

- **Parallel Execution**: The algorithm is parallelized using Pthreads to leverage multi-core processors for faster contour extraction.
  
- **Efficient Resource Management**: Dynamic memory allocation and deallocation are optimized to prevent memory leaks and ensure efficient memory usage during execution.

- **Scalability**: The implementation is designed to scale efficiently with the number of available processor cores, enabling performance improvements on multi-core systems.

## Usage

To use the parallelized Marching Squares algorithm:

1. Compile the source code using a C compiler with Pthreads support.
2. Execute the compiled binary, providing the input image file path and the desired output file path as command-line arguments.

Example:

```
./marching_squares input_image.ppm output_contour.ppm
```
