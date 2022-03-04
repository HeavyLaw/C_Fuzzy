#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>


#define MAT_LEGAL_CHECKING
#define MAT_INIT_FLAG 0x5C

typedef struct{
	int row,col;
	float **element;
	unsigned char init;
}Mat;

Mat* MatCreate(Mat* mat,int row,int col);
void MatDelete(Mat* mat);
Mat* MatSetVal(Mat* mat,float* val);
void MatDump(const Mat* mat);
float FuzzyIst(float a,float b);
float FuzzyCbn(float a,float b);

Mat* MatCreate(Mat* mat,int row,int col)
{
	int i;

#ifdef MAT_LEGAL_CHECKING
	if(mat->init == MAT_INIT_FLAG){
		if(mat->row == row && mat->col == col)
			return mat;
		else
			MatDelete(mat);
	}
#endif

	mat->element = (float**)malloc(row * sizeof(float*));
	for(i = 0;i<row;i++){
		mat->element[i] = (float*)malloc(col * sizeof(float));
	}

	if(mat->element == NULL){
		return NULL;
	}
	mat->row = row;
	mat->col = col;
	mat->init = MAT_INIT_FLAG;
	
	return mat;
}

void MatDelete(Mat* mat)
{
	int i;

#ifdef MAT_LEGAL_CHECKING
	if(mat->init == MAT_INIT_FLAG){
		return ;
	}
#endif
	
	for(i = 0;i < mat->row; i++)
		free(mat->element[i]);
	free(mat->element);

	mat->init = 0;
}

Mat* MatSetVal(Mat* mat,float* val)
{
	int row,col;

	for(row = 0;row < mat->row; row++){
		for(col = 0;col < mat->col; col++){
			mat->element[row][col] = val[col + row * mat->col];
		}
	}

	return mat;
}

void MatDump(const Mat* mat)
{
	int row,col;
#ifdef MAT_LEGAL_CHECKING
	if(mat == NULL){
			return ;
	}
#endif

	printf("Mat %dx%d:\n",mat->row,mat->col);
	for(row = 0;row < mat->row; row++){
		for(col = 0;col < mat->col; col++){
			printf("%.2f\t",mat->element[row][col]);
		}
		printf("\n");
	}
}

Mat* MatTrans(Mat* src, Mat* dst)
{
    int row, col;

#ifdef MAT_LEGAL_CHECKING
    if( src->row != dst->col || src->col != dst->row ){
        printf("err check, unmatch matrix for MatTranspose!\n");
        MatDump(src);
        MatDump(dst);
        return NULL;
    }
#endif

    for(row = 0 ; row < dst->row ; row++){
        for(col = 0 ; col < dst->col ; col++){
            dst->element[row][col] = src->element[col][row];
        }
    }

    return dst;
}

bool MatEqual(Mat* mat1,Mat* mat2){
	int row,col;
	bool flag;

#ifdef MAT_LEGAL_CHECKING
	if(mat1->row != mat2->row || mat1->col != mat2->col){
		printf("err check,unmatch matrix for MatEqual!\n");
		MatDump(mat1);
		MatDump(mat2);
		return 0;
	}
#endif

	for(row = 0;row < mat1->row;row++){
		for(col = 0;col < mat1->col;col++){
			if(mat1->element[row][col] == mat2->element[row][col])
				flag = 1;
			else
			{
				flag = 0;
				goto here;
			}
		}
	}
	here:
	return flag;
}

void MatSwap(Mat *mat1,Mat *mat2){
	int row,col;
	Mat mat3;

#ifdef MAT_LEGAL_CHECKING
	if(mat1->row != mat2->row || mat1->col != mat2->col){
		printf("err check,unmatch matrix for MatSwap!\n");
		MatDump(mat1);
		MatDump(mat2);
	}
#endif

	MatCreate(&mat3,mat1->row,mat1->col);
	for(row = 0;row < mat1->row;row++){
		for(col = 0;col <mat1->col;col++){
			mat3.element[row][col] = mat1->element[row][col];
			mat1->element[row][col] = mat2->element[row][col];
			mat2->element[row][col] = mat3.element[row][col];
		}
	} 
}

Mat* FuzzyMatMul(Mat* src1,Mat* src2,Mat* dst){
	 int row,col;
	 int i;
	 float temp;

#ifdef MAT_LEGAL_CHECKING
	 if(src1->col != src2->row || src1->row != dst->row || src2->col != dst->col){
		 printf("err check,unmatch matrix for FuzzyMatMul!\n");
		 MatDump(src1);
		 MatDump(src2);
		 MatDump(dst);
		 return NULL;
	 }
#endif

	 for(row = 0; row < dst->row; row++){
		 for(col = 0; col < dst->col; col++){
			 temp = 0.0f;
			 for(i = 0; i < src1->col; i++){
				 temp = FuzzyCbn(temp,FuzzyIst(src1->element[row][i],src2->element[i][col]));
			 }
			 dst->element[row][col] = temp;
			 printf("%.2f\t",dst->element[row][col]);
		 }
		 printf("\n");
	 }
	 return dst;
}

Mat* MaxNormal(Mat* dst,Mat* src)
{
	int row,col;
	float *max = (float*)malloc(src->col * sizeof(float));

#ifdef MAT_LEGAL_CHECKING
	if(src == NULL){
			return NULL;
	}
#endif
	
	printf("\n");
	for(col = 0;col < src->col; col++){
		max[col] = src->element[0][col];
		for(row = 1;row < src->row; row++){
			if(src->element[row][col] > max[col])
				max[col] = src->element[row][col];
			else
				max[col] = max[col];
		}
		printf("Maximum value of the %dth column:%.0f\n",col+1,max[col]);
	}

	printf("\nMaximum normalized matrix X:\n");	
	for(row = 0;row < dst->row; row++){
		for(col = 0;col < dst->col; col++){
			dst->element[row][col] = src->element[row][col]/max[col];
			printf("%.2f\t",dst->element[row][col]);
		}
		printf("\n");
	}

	return dst;
}

float FuzzyIst(float a,float b){
	if(a >= b)
		a = b;
	else
		a = a;

	return a;
}

float FuzzyCbn(float a,float b){
	if(a >= b)
		a = a;
	else
		a = b;

	return a;
}

Mat* MaxMin(Mat* dst,Mat* src,Mat* src_T){
	int row,col;
	int k;	
	float temp1,temp2;

#ifdef MAT_LEGAL_CHECKING
	if(src == NULL){
			return NULL;
	}
#endif

	printf("\nConstruction of Fuzzy Similarity Matrix by Maximum and Minimum Method R:\n");	
	for(row=0;row < dst->row;row++){
		for(col=0;col < dst->col;col++){
			temp1 = 0;
			temp2 = 0;
			for(k=0;k < src->col;k++){
				temp1 += FuzzyIst(src->element[row][k],src_T->element[k][col]);
				temp2 += FuzzyCbn(src->element[row][k],src_T->element[k][col]);
			}
			dst->element[row][col] = temp1/temp2;
			printf("%.2f\t",dst->element[row][col]);
		}
		printf("\n");
	}	

	return dst;
}

int FuzzyLambda(Mat* mat,Mat* lambda){
	int row,col;
	float *temp = (float*)malloc((mat->col * (mat->row-1)/2)* sizeof(float));
	int i,j,k = 0,lambda_i = 1;
	float tmp;

#ifdef MAT_LEGAL_CHECKING
	for(int i2 = 0;i2 < mat->row;i2++){
		if(mat == NULL || mat->element[i2][i2] != 1){
			printf("err check,unmatch matrix for printLambda!\n");
			return NULL;
		}
	}
#endif
	
	printf("\nArrange the elements of t(R) from large to small:\n");
	for(row = 0;row < mat->row-1;row++){
		for(col = row + 1;col < mat->col;col++){
				temp[k++] = mat->element[row][col];
		}
	}
	
	for(i = 1;i < mat->col * (mat->row-1)/2;i++){
		for(j = 0;j < i; j++){
			if(temp[j] == temp[i])
			{
				temp[i] = 100.00f;
				break;
			}
		}
	}
	for(i = 1;i < mat->col * (mat->row-1)/2;i++){
		for(j = 0;j < i; j++){
			if(temp[j] < temp[i])
			{
				tmp = temp[i];
				temp[i] = temp[j];
				temp[j] = tmp;
			}
		}
	}
	lambda->element[0][0] = 1;
	printf("%.2f\t",lambda->element[0][0]);
	for(k = 0;k < mat->col * (mat->row-1)/2;k++){
		if(temp[k] != 100.00f)
		{
			printf("%.2f\t",temp[k]);
			lambda->element[0][lambda_i++] = temp[k];
		}
	}

	printf("\n");
	return lambda_i;
}

void FuzzyCluster(float lambda,Mat *mat){
	int row,col;
	Mat temp;
	int i,flag;

#ifdef MAT_LEGAL_CHECKING
	if(mat == NULL){
		return ;
	}
#endif
	MatCreate(&temp,mat->row,mat->col);
	
	printf("Take %.2f=1 to get T (%.2f) :\n",lambda,lambda);
	for(row = 0;row < mat->row;row++){
		for(col = 0;col < mat->col;col++){
			if(mat->element[row][col] >= lambda)
			{
				temp.element[row][col] = 1;
				printf("%.0f\t",temp.element[row][col]);
			}
			else
			{
				temp.element[row][col] = 0;
				printf("%.0f\t",temp.element[row][col]);
			}
		}
		printf("\n");
	}
	for(row = 0;row < temp.row;row++){
		flag = 0;
		for(i = 0;i < row;i++){
			for(col = 0;col < temp.col;col++){
			if(temp.element[row][col] != temp.element[i][col])
				break;
			}
			if(col == temp.col)
			{
				printf("x%d\t",i+1);
				flag++;
			}			
		}
		if(flag != 0)
			printf("\nx%d\n",row+1);
	}

}	


void main(){
	Mat a;
	Mat X,X_T;
	Mat R;
	Mat T;
	Mat Lambda;
	int n,i = 1;

	float val[] = {80,10,6,2,
					50,1,6,4,
					90,6,4,6,
					40,5,7,3,
					10,1,2,4};

	MatCreate(&a,5,4);
	MatCreate(&X,5,4);
	MatCreate(&R,5,5);
	MatCreate(&T,5,5);
	MatCreate(&X_T,X.col,X.row);
	MatCreate(&Lambda,1,T.col * (T.row-1)/2);

	MatSetVal(&a,val);
	MatDump(&a);

	MaxNormal(&X,&a);
	MatTrans(&X,&X_T);
	MaxMin(&R,&X,&X_T);

	while(!MatEqual(&R,&T)){
		printf("Transitive Closure with Square Method R(%d):\n",++i);
		FuzzyMatMul(&R,&R,&T);
		MatSwap(&R,&T);
	}

	n = FuzzyLambda(&T,&Lambda);
	for(int j = 0;j < n;j++)
	{
		FuzzyCluster(Lambda.element[0][j],&T);
	}
}



