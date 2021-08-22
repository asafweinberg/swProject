#include "spkmeans.c"

void printMat(matrix *m);

int main(int argc,char* argv[]){

    matrix* m1= newMatrix(3,3);
    matrix* m2= newMatrix(3,3);
    matrix* m3, *m5;
    // m1->data[1][1]=2.5;
    // m1->data[1][2]=3.5;
    // m1->data[1][3]=4.5;
    // m1->data[2][1]=5.5;
    // m1->data[2][2]=-6.5;
    // m1->data[2][3]=7.5;
    // m1->data[3][1]=-28;
    // m1->data[3][2]=1;
    // m1->data[3][3]=0;

    // m2->data[1][1]=-6.9;
    // m2->data[1][2]=4.87;
    // m2->data[1][3]=9;
    // m2->data[2][1]=1;
    // m2->data[2][2]=2.2;
    // m2->data[2][3]=3.3;
    // m2->data[3][1]=6.6;
    // m2->data[3][2]=7.7;
    // m2->data[3][3]=8.8;
 //l   printf("1");
    m1->data[0][0]=1;
    m1->data[0][1]=2;
    m1->data[0][2]=3;
    m1->data[1][0]=4;
    m1->data[1][1]=5;
    m1->data[1][2]=6;
    m1->data[2][0]=7;
    m1->data[2][1]=8;
    m1->data[2][2]=9;
    
    

    m2->data[1][1]=1;
    m2->data[1][2]=1;
    m2->data[1][0]=1;
    m2->data[2][1]=1;
    m2->data[2][2]=1;
    m2->data[2][0]=1;
    m2->data[0][1]=1;
    m2->data[0][2]=1;
    m2->data[0][0]=1;

    m3=formMatI(3);
    // printMat(m1);
    // printMat(m3);
 //   printf("2");

    //matrix * m3=addMatrices(m1,m2,true,true,true);
    // matrix * m4=mulMatrices(m1,m3,true,true);
  //  printf("3");

    // printMat(m4);
    // free(m1);
    // printf("%p",m1);
    // printf("4gg");
    m1=minusRootMat(m1,true);
    printMat(m1);
}


