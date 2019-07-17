#include <cstdio> // printf
// from GROMACS
void calc_triclinic_images(double (*img)[4], double const box[3][4]) {
  
    int constexpr XX=0, YY=1, ZZ=2;
    /* Calculate 3 adjacent images in the xy-plane */
    for(int d = 0; d < 3; ++d) img[0][d] = box[XX][d];
    for(int d = 0; d < 3; ++d) img[1][d] = box[YY][d];
    if (img[1][XX] < 0) {
        for(int d = 0; d < 3; ++d) img[1][d] *= -1;
    }
    for(int d = 0; d < 3; ++d) img[2][d] = img[1][d] - img[0][d];

    /* Get the next 3 in the xy-plane as mirror images */
    for (int i = 0; i < 3; ++i) {
        for(int d = 0; d < 3; ++d) img[3+i][d] = -img[i][d];
    } // i

    /* Calculate the first 4 out of xy-plane images */
    for(int d = 0; d < 3; ++d) img[6][d] = box[ZZ][d];
    if (img[6][XX] < 0) {
        for(int d = 0; d < 3; ++d) img[6][d] *= -1;
    }
    for (int i = 0; i < 3; ++i) {
        for(int d = 0; d < 3; ++d) img[7+i][d] = img[6][d] + img[i+1][d];
    } // i

    /* Mirror the last 4 from the previous in opposite rotation */
    for (int i = 0; i < 4; ++i) {
        for(int d = 0; d < 3; ++d) img[10+i][d] = -img[6 + (2+i) % 4][d];
    } // i
}

int main(){
  
  double img[16][4], box[3][4] = {{1,0,0,  0}, {0,2,0,  0}, {0,0,3,  0}};
  calc_triclinic_images(img, box);
  for(int i = 0; i < 14; ++i) {
      printf("# image %2d     %.6f  %.6f  %.6f\n", i, img[i][0], img[i][1], img[i][2]);
  }
  return 0;
} // main
