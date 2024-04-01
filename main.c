/*************************************************
*************CUBIC SPLINE PROGRAM*****************
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
//
#define MAX_LEN_BUFF 100
#define MAX_POINTS 50
#define MAX_PLOTS 5
#define MAX_ITERATIONS 5000
#define H 0.05
#define EPS 1e-7
#define DIR_MODE 0777
#define MODE_MAX 1
#define MODE_MIN -1
//
#define GET_VAR_NAME(var)(#var)

const char *PLOT_COMMANDS_PATH = "plt_commands.txt";
const char *COLORS[] = {"orange", "blue", "green", "gold", "purple", "black"};

typedef struct Point {
    double x;
    double y;
} Point;

typedef struct ComplexFunc {
    Point *v1, *v2;
    int size_1, size_2;
    double **coef_1, **coef_2;
    double left, right;
} ComplexFunc;

void gen_spline(const Point *vec_p, const int v_size, double **coef_mat);
double interpolate(const Point *vec_p, const int v_size, double **coef_mat, double desired_x);
void find_intersection_of_sp(ComplexFunc cmf);
double find_min_distance_btw_sp(ComplexFunc cmf);
void printAnswer(int flag, ComplexFunc cmf);

double pow(double base, double power)
{
    double res = base;
    for(int i = 0; i < power - 1; i++)
        res *= base;
    return res;
}

double fabs(double x)
{
    return x >= 0 ? x : -x;
}

double sqrt(double input)
{
    const double SQRT_EPS = 1e-15;
    double nx;
    double x = 1;
    for(;;) {
	    nx = (x + input / x) / 2;
	    if (fabs(x - nx) < SQRT_EPS)  
            break;
	    x = nx;
    }
    return x;
}

double compare(double a, double b, int mode)
{
    if(mode == MODE_MAX)
        return a > b ? a : b;  // max
    return a < b ? a : b;      // min
}

void init_vec(Point *vec, const int size, const char *name)
{
    printf("Init points for the \"%s\" plot separated by space (x,y): \n", name);
    for(int i = 0; i < size; i++)
        scanf("%lf %lf", &vec[i].x, &vec[i].y);
}

void init_paths(char **paths, const int size)
{
    int ind;
    mkdir("data/", DIR_MODE);
    for(int i = 0; i < size; i += 2)
    {
        ind = i/2 + 1;
        char buff[MAX_LEN_BUFF];

        sprintf(buff, "data/func_%d/", ind);
        mkdir(buff, DIR_MODE);

        sprintf(paths[i], "data/func_%d/base_points_%d.txt", ind, ind);
        sprintf(paths[i + 1], "data/func_%d/line_%d.txt", ind, ind);
    } 
}

void init_data_func_file(Point *vec_p, const int v_size, char *path_vert, char *path_line, double **coef)
{
    FILE *file_vert = fopen(path_vert, "w+");
    FILE *file_line = fopen(path_line, "w+");

    // data-vertexes
    for(int i = 0; i < v_size; i++)
    {
        char buff[MAX_LEN_BUFF];
        sprintf(buff, "%lf  %lf\n", vec_p[i].x, vec_p[i].y);
        fputs(buff, file_vert);
    }
    fclose(file_vert);

    // data-line
    gen_spline(vec_p, v_size, coef);
    for(double x = vec_p[0].x; x <= vec_p[v_size - 1].x; x += H)
    {
        char buff[MAX_LEN_BUFF];
        sprintf(buff, "%lf  %0.9lf\n", x, interpolate(vec_p, v_size, coef, x));
        fputs(buff, file_line);
    }
    fclose(file_line);
}

void display_vec(const Point *vec, const int size, const char *name)
{
    printf("\n-> vector \"%s\":\n", name);
    for(int i = 0; i < size; i++)
        printf("( %lf ; %lf )\n", vec[i].x, vec[i].y);
    printf("-----\n");
}

void display_plots(char **paths, int size)
{

    int ind;
    char buff[MAX_LEN_BUFF];
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    FILE *pre_instruction = fopen(PLOT_COMMANDS_PATH, "r");

    while(fgets(buff, MAX_LEN_BUFF, pre_instruction) != NULL)
        fprintf(gnuplotPipe, "%s\n", buff); 

    char *cmd_plot = malloc(sizeof(char)*MAX_LEN_BUFF*size);
    strcat(cmd_plot, "plot");

    for(int i = 0; i < size; i += 2)
    {
        ind = i/2 + 1;
        char plot_buff[MAX_LEN_BUFF], cmd_style[MAX_LEN_BUFF];
        sprintf(cmd_style, "set style line %d lc rgb \"%s\" lw 3 pt 6", ind, COLORS[ind - 1]);
        sprintf(plot_buff, " \"%s\" w l title \"f_%d\" ls %d, \"%s\" w p notitle ls %d", 
                paths[i + 1], ind, ind, paths[i], ind);   // -> w p/l/pl
        
        if(i + 2 < size)
            strcat(plot_buff, ",");
        strcat(cmd_plot, plot_buff);
        fprintf(gnuplotPipe, "%s\n", cmd_style);
    }
    fprintf(gnuplotPipe, "%s\n", cmd_plot);
    
    fflush(gnuplotPipe);
    fclose(gnuplotPipe);
    free(cmd_plot);
}

/*
    Performs Gauss-Elimination and returns the Upper triangular matrix and solution of equations
*/
void gauss_elimination(const int n_iter, double *sig_tmp, double **tdm)
{
    for(int i = 0; i < n_iter - 2; i++)
    {
        // begin Gauss elimination
        for(int k = i + 1; k < n_iter - 1; k++)
        {
            double term = tdm[k][i]/tdm[i][i];
            for(int j = 0; j < n_iter; j++)
                tdm[k][j] = tdm[k][j] - term*tdm[i][j];
        }
    }
    // begin back-substitution
    for(int i = n_iter - 2; i >= 0; i--)
    {
        sig_tmp[i] = tdm[i][n_iter - 1];
        for(int j = i + 1; j< n_iter - 1; j++)
        {
            sig_tmp[i] = sig_tmp[i] - tdm[i][j]*sig_tmp[j];
        }
        sig_tmp[i] = sig_tmp[i]/tdm[i][i];
    }      
}

/*
    Cubic spline coefficients calculator
*/
void calc_coef_sp(const Point *p_vec, const int num_iter, double *h, double *sig, double **coef)
{
    for(int i = 0; i < num_iter; i++)
    {
        coef[i][3] = p_vec[i].y;
        coef[i][1] = sig[i]/2.0;
        coef[i][0] = (sig[i + 1] - sig[i])/(h[i]*6.0);
        coef[i][2] = (p_vec[i + 1].y - p_vec[i].y)/h[i] - h[i]*(2*sig[i] + sig[i + 1])/6.0;
    }
}

/*
    Generate the tridiagonal augmented matrix for cubic spline for equidistant data-points
*/
void gen_tgm(const Point *p_vec, const int num_iter, double *h, double **tdm)
{
    for(int i = 0; i < num_iter - 1; i++)
        tdm[i][i] = 2*(h[i] + h[i + 1]);

    for(int i = 0; i < num_iter - 2; i++)
    {
        tdm[i][i + 1] = h[i + 1];
        tdm[i + 1][i] = h[i + 1];
    }
    for(int i = 1; i < num_iter; i++)
        tdm[i - 1][num_iter - 1] = (p_vec[i + 1].y - p_vec[i].y)*6/(double)h[i] - 
                                   (p_vec[i].y - p_vec[i - 1].y)*6/(double)h[i - 1];
} 

double interpolate(const Point *vec_p, const int v_size, double **mat_coef, double desired_x)
{
	int i = 0;
    if(desired_x > vec_p[v_size - 1].x || desired_x < vec_p[0].x)
    {
        printf("The point is outside the interval! (%lf)\n\n", desired_x);
        exit(0);
    }
	while(vec_p[i].x < desired_x)
		i++;
	i--;
    if(i < 0) i = 0;
	return  mat_coef[i][0]*pow(desired_x - vec_p[i].x, 3) + // A
            mat_coef[i][1]*pow(desired_x - vec_p[i].x, 2) + // B
            mat_coef[i][2]*(desired_x - vec_p[i].x) +       // C
            mat_coef[i][3];                                 // D
}

void gen_spline(const Point *vec_p, const int v_size, double **mat_coef)
{
    double h[v_size - 1];
    double sig[v_size]; 
    double sig_tmp[v_size - 2];

    for(int i = 0; i < v_size - 1; i++)
        h[i] = vec_p[i+1].x - vec_p[i].x;

    sig[0] = 0;
    sig[v_size - 1] = 0;

	double **tdm = malloc(sizeof(double*)*(v_size - 2));
	for(int i = 0; i < v_size - 2; i++)
		tdm[i] = malloc(sizeof(double)*(v_size - 1)); 
    
    gen_tgm(vec_p, v_size - 1, h, tdm);  
    gauss_elimination(v_size - 1, sig_tmp, tdm);

    for(int i = 1; i < v_size - 1; i++)
        sig[i] = sig_tmp[i-1];
    
    calc_coef_sp(vec_p, v_size - 1, h, sig, mat_coef);
    free(tdm);
}

void free_space(const int n, ...)
{
    va_list args;
    va_start(args, n);
    for(int i = 0; i < n; i++)
    {
        void *va = va_arg(args, void*);
        free(va);
    }
    va_end(args);
}

int is_natural(const int num)
{
    if(num < 1)
    {
        printf("Set a natural number!\n\n");
        exit(0);
    }
    return 1;
}

int main()
{
    int num_plots;
    printf("Define num of plots: ");
    scanf("%d", &num_plots);

    if(is_natural(num_plots) && num_plots > 5)
    {
        printf("Max size of plots is %d!\n\n", MAX_PLOTS);
        return 0;
    }

    double ***coef = malloc(sizeof(double**)*num_plots);

    char **paths = malloc(sizeof(char*)*(2*num_plots + 1));
    for(int i = 0; i <= 2*num_plots; i++)
        paths[i] = i != 2*num_plots ? malloc(sizeof(char)*MAX_LEN_BUFF) : NULL;

    init_paths(paths, 2*num_plots);

    int dim_plots[num_plots];
    printf("Set the dimension of the vectors for plots separated by space: ");
    for(int i = 0; i < num_plots; i++)
    {
        scanf("%d", (dim_plots + i));
        is_natural(dim_plots[i]);
        coef[i] = malloc(sizeof(double*)*(dim_plots[i] - 1));
        for(int j = 0; j < dim_plots[i]; j++)
            coef[i][j] = malloc(sizeof(double)*4);
    }
        
    Point **arr_vp = malloc(sizeof(Point*)*num_plots);
    for(int i = 0; i < num_plots; i++)
    {
        char buff[MAX_LEN_BUFF];
        arr_vp[i] = malloc(sizeof(Point)*dim_plots[i]);
        sprintf(buff, "%d-vec", i + 1);

        init_vec(arr_vp[i], dim_plots[i], buff);
        init_data_func_file(arr_vp[i], dim_plots[i], paths[2*i], paths[2*i + 1], coef[i]);
    }
    display_plots(paths, 2*num_plots);


    if(num_plots > 1)
    {
        int ind_1, ind_2;
        printf("\nEnter two plot-index to find some intersection (if exists), separated by space (from 1 to %d): ", num_plots);
        scanf("%d %d", &ind_1, &ind_2);
        ind_1--; ind_2--;

        ComplexFunc cmf = { .v1 = arr_vp[ind_1], .v2 = arr_vp[ind_2],
                              .size_1 = dim_plots[ind_1], .size_2 = dim_plots[ind_2],
                              .coef_1 = coef[ind_1], .coef_2 = coef[ind_2],
                              .left = compare(arr_vp[ind_1][0].x, arr_vp[ind_2][0].x, MODE_MAX),
                              .right = compare(arr_vp[ind_1][dim_plots[ind_1] - 1].x, 
                                               arr_vp[ind_2][dim_plots[ind_2] - 1].x, MODE_MIN) };
        
        if(arr_vp[ind_1][0].x > arr_vp[ind_2][dim_plots[ind_2] - 1].x || arr_vp[ind_2][0].x > arr_vp[ind_1][dim_plots[ind_1] - 1].x)
            printAnswer(1, cmf);
        else
            find_intersection_of_sp(cmf);
    }

    free_space(3, paths, coef, arr_vp);
    return 0;
}

/*************************************************
*************SPLINES INTERSECTION(S)**************
*************************************************/

double *x, *x_next, *fx, **mat_Jacobi;

void initFields(int n, double l, double r)
{
    x = malloc(n*sizeof(double));
    x_next = malloc(n*sizeof(double));
    fx = malloc(n*sizeof(double));
    mat_Jacobi = malloc(n*sizeof(double*));
    for(int i = 0; i < n; i++)
    {
        mat_Jacobi[i] = malloc(sizeof(double)*n);
        x[i] = (l+r)/2;
    }
}

double get_func_val(int pos, ComplexFunc cmf)
{
    if(x[0] > cmf.right)
        x[0] = cmf.right - EPS;
    if(x[0] < cmf.left)
        x[0] = cmf.left + EPS;

    if(pos) 
        return interpolate(cmf.v2, cmf.size_2, cmf.coef_2, x[0]) - x[1];
    return interpolate(cmf.v1, cmf.size_1, cmf.coef_1, x[0]) - x[1];
}

void getCofactor(double **mat, double **temp, int p, int q, int size)
{
    int i = 0, j = 0;
    for(int row = 0; row < size; row++)
    {
        for(int col = 0; col < size; col++)
        {
            if(row != p && col != q)
            {
                temp[i][j] = mat[row][col];
                j++;
                if(j == size - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double determinant(double **mat, int n)
{
    int i, j, k, factor = 1;
    double det = 0, **new_mat;
    if(mat == NULL)
    {
        printf("\nError in det calculation!\n");
        exit(0);
    }
    if(n == 1)
        return mat[0][0];
    for(i = 0; i < n; i++)
    {
        if(NULL == (new_mat = malloc((n - 1)*sizeof(double*))))
            return -1;
        for(j = 0; j < n - 1; j++)
            if(NULL == (new_mat[j] = malloc((n - 1)*sizeof(double))))
                return -1;

        for(j = 1; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                if(k == i) continue;
                new_mat[j - 1][k < i ? k : (k - 1)] = mat[j][k];
            }
        }
        det += factor * mat[0][i] * determinant(new_mat, n - 1);
        factor *= -1;
        free(new_mat);
    }
    return det;
}

void adjointMat(double **mat, double **adj, int size)
{
    if(size == 1)
    {
        adj[0][0] = 1;
        return;
    }
    int sign = 1;
    double **temp = malloc(size*sizeof(double*));
    for(int i = 0; i < size; i++)
        temp[i] = malloc(size*sizeof(double));;

    for(int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            getCofactor(mat, temp, i, j, size);
            sign = ((i + j)%2 == 0) ? 1 : -1;
            adj[j][i] = sign*determinant(temp, size - 1);
        }
    }
    free(temp);
}

int is_identical(ComplexFunc cmf)
{
    for(int i = 0; i < compare(cmf.size_1, cmf.size_2, MODE_MIN); i++)
    {
        if(cmf.v1[i].x != cmf.v2[i].x)
            return 0;
    }
    return 1;
}

double** inverseMat(double **mat, int size, ComplexFunc cmf)
{
    double det = determinant(mat, size);
    if(det == 0)
    {
        if(is_identical(cmf))
        {
            printf("\nThese splines have an inf number of intersection points!\n\n");
            exit(0);
        }
        else
            printAnswer(1, cmf);
    }
    double **adj = malloc(size*sizeof(double*)),
           **inv = malloc(size*sizeof(double*));

    for(int i = 0; i < size; i++)
    {
        adj[i] = malloc(size*sizeof(double));
        inv[i] = malloc(size*sizeof(double));
    }
    adjointMat(mat, adj, size);

    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            inv[i][j] = adj[i][j] / det;

    free(adj);
    return inv;
}

double getPartialDerivative(int pos_expr, int pos_x, int size, ComplexFunc cmf)
{
    *(x + pos_x) += H;    
    const double h1 = get_func_val(pos_expr, cmf);

    *(x + pos_x) -= 2*H;    
    const double h2 = get_func_val(pos_expr, cmf);

    *(x + pos_x) += H;
    return (h1 - h2)/(2*H);
}

void fillMatJacobi(double **mat_Jacobi, int size, ComplexFunc cmf)
{
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            mat_Jacobi[i][j] = getPartialDerivative(i, j, size, cmf);
}

int iterCriterion(int size)
{
    for(int i = 0; i < size; i++)
        if((fx[i] > 0 ? fx[i] : -fx[i]) > EPS)
            return 1;
    return 0;
}

double* diffVec(double* a, double* b, int size)
{
    for(int i = 0; i < size; i++)
        *(a + i) -= *(b + i);
    return a;
}

void fillEquationVec(int size, ComplexFunc cmf)
{
    for(int i = 0; i < size; i++)
        *(fx + i) = get_func_val(i, cmf);
}

double* getMatVecMultiplication(double **mat, double *vec, int size)
{
    double sm, *d_vec;
    d_vec = malloc(size*sizeof(double));

    for(int i = 0; i < size; i++)
    {
        sm = 0;
        for(int j = 0; j < size; j++)
            sm += mat[i][j]*vec[j];
        *(d_vec + i) = sm;
    }
    return d_vec;
}

void nextStep(int n, ComplexFunc cmf)
{
    fillMatJacobi(mat_Jacobi, n, cmf);
    fillEquationVec(n, cmf);
    x_next = diffVec(x, getMatVecMultiplication(inverseMat(mat_Jacobi, n, cmf), fx, n), n);
}

void printAnswer(int flag, ComplexFunc cmf)
{
    if(!flag)
        printf("\nThe point of intersection is: (%lf ; %lf)\n\n", x[0], x[1]);
    else
    {
        printf("\nThese splines do not intersect!\n");
        printf("Min distance is %lf\n\n", find_min_distance_btw_sp(cmf));
    }
    printf("---------\n");
}

void find_intersection_of_sp(ComplexFunc cmf)
{
    int n = 2, steps = 0;
    initFields(n, cmf.left, cmf.right);
    nextStep(n, cmf);

    while(iterCriterion(n) && steps <= MAX_ITERATIONS)
    {
        x = x_next;
        nextStep(n, cmf);
        steps++;
    }
    printAnswer((steps > MAX_ITERATIONS && iterCriterion(n)), cmf);
    free_space(2, x, fx);
    free(mat_Jacobi);
}

/*************************************************
*************MIN DISTANCE B/W SPLINES*************
*************************************************/
 
#define COEF 1.1

typedef struct
{
   double x1, x2;
} vector;
 
double obj_func(vector v, ComplexFunc cmf)
{
    double y1, y2;

    if(v.x1 > cmf.v1[cmf.size_1 - 1].x)
        v.x1 = cmf.v1[cmf.size_1 - 1].x - EPS;
    if(v.x1 < cmf.v1[0].x)
        v.x1 = cmf.v1[0].x + EPS;
    
    if(v.x2 > cmf.v2[cmf.size_2 - 1].x)
        v.x2 = cmf.v2[cmf.size_2 - 1].x - EPS;
    if(v.x2 < cmf.v2[0].x)
        v.x2 = cmf.v2[0].x + EPS;

    y1 = interpolate(cmf.v1, cmf.size_1, cmf.coef_1, v.x1);
    y2 = interpolate(cmf.v2, cmf.size_2, cmf.coef_2, v.x2);
    return sqrt(pow(v.x2 - v.x1, 2) + pow(y2 - y1, 2));
}
 
vector gradient(vector v)
{
   vector grad;
   grad.x1 = COEF*v.x1;
   grad.x2 = COEF*v.x2;
   return grad;
}
 
double make_simple_fx(double x, vector grad, vector xj, ComplexFunc cmf)
{
   vector buffer;
 
   buffer.x1 = xj.x1 - x*grad.x1;
   buffer.x2 = xj.x2 - x*grad.x2;
 
   return obj_func(buffer, cmf);
}
 
double golden_selection(double a, double b, vector gradient, vector x, ComplexFunc cmf)
{
    const double fi = 1.6180339887;
    double x1, x2, y1, y2;
 
    x1 = b - ((b-a)/fi);
    x2 = a + ((b-a)/fi);
    y1 = make_simple_fx(x1, gradient, x, cmf);
    y2 = make_simple_fx(x2, gradient, x, cmf);

    while(fabs(b - a) > EPS)
    {
        if(y1 <= y2)
        {
            b = x2;
            x2 = x1;
            x1 = b - (b - a)/fi;
            y2 = y1;
            y1 = make_simple_fx(x1, gradient, x, cmf);
        }
        else
        {
            a = x1;
            x1 = x2;
            x2 = a + (b - a)/fi;
            y1 = y2;
            y2 = make_simple_fx(x2, gradient, x, cmf);
        }
    }
    return (a + b)/2;
}
 
vector calculate(vector x, vector gradient, double lambda)
{
   vector buffer;
 
   buffer.x1 = x.x1 - lambda*gradient.x1;
   buffer.x2 = x.x2 - lambda*gradient.x2;
 
   return buffer;
}
 
vector grad_down(vector x, ComplexFunc cmf)
{
    vector current = x;
    vector last;
    int steps = 0;
    do
    {
        last = current; 
        vector grad = gradient(current); 
        double lambda = golden_selection(0, 0.05, grad, current, cmf); 
        current = calculate(current, grad, lambda); 
        steps++;

    } while(fabs(obj_func(current, cmf) - obj_func(last, cmf)) > EPS && steps <= MAX_ITERATIONS*50);
    
    return current; 
}

double find_min_distance_btw_sp(ComplexFunc cmf)
{
    vector init_approximation = {.x1 = (cmf.v1[0].x + cmf.v1[cmf.size_1 - 1].x)/2,
                                 .x2 = (cmf.v2[0].x + cmf.v2[cmf.size_2 - 1].x)/2 };
    vector rez = grad_down(init_approximation, cmf);

    return obj_func(rez, cmf);
}