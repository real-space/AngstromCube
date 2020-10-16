#pragma once

template<typename real_t>
int dump_to_file(
    char const *filename, 
    int const N, // number of rows
    real_t const y_data[], // 2D-data array [N*Stride]
    double *x_axis=nullptr, // x-axis[N], replaced by 0,1,2,... if 0x0
    int const Stride=1, // Stride, default=1 for 1D functions
    int const M=1, // number of columns to display
    char const *title=nullptr, // title line in the file
    int const echo=0) // report writing to stdout, 0: suppress report
{
    FILE *f = fopen(filename, "w");
    if (nullptr == f) {
        if (echo > 1) printf("# %s Error opening file %s!\n", __func__, filename);
        return 1;
    } // failed to open

    fprintf(f, "# %s\n", title); // print this line also if title==nullptr

    for(int i = 0; i < N; i++) {
        if (nullptr != x_axis) {
            fprintf(f, "%g ", x_axis[i]);
        } else {
            fprintf(f, "%d ", i);
        } // x_axis given
        for(int j = 0; j < M; ++j) {
            fprintf(f, "%g ", y_data[i*Stride + j]);
        } // j
        fprintf(f, "\n"); // new line
    } // i
    fprintf(f, "\n"); // new line at the end of the file

    fclose(f);
    if (echo > 3) printf("# file %s written with %d x %d (of %d) data entries.\n", filename, N, M, Stride);
    return 0;
} // dump_to_file

#define here if (echo > 5) { printf("\n# here: %s %s:%i\n\n", __func__, __FILE__, __LINE__); fflush(stdout); }
