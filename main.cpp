#include "NTL_head.h"
#include "gmp_head.h"
#include <algorithm>
#include <sys/types.h>
#include <unistd.h>

#define DIM 125
#define EQU 200
#define LEN 256
#define BITS 224
#define UNKNOWN 222
#define COEFF "Coeff_8.txt"
#define KNOWN "Known_8.txt"
#define BKZ_BLOCK_SIZE 70
#define BKZ_SET_MAX_TIME 12000
#define TEST_TIMES 1000000000
#define PRINT_HIDE 0
#define PRINT_HADAMARD 1
#define PARALLEL_TIMES 8
#define PROCEEDS 47
double noise[DIM + 1];
int computed = 0;
double noise_times = 0.2;

struct TestData {
	int data_id;
	double reduce_time;
	int success_times;
	int all_times;
	double average_nv_time;
  RR hadamard;
};
ostream& operator << (ostream& out, TestData& td)
{
	out << td.data_id << '\t' << td.reduce_time << '\t';
 /*
	if (td.success_times == 0)
		out << "<1";
	else
		out << td.success_times;
	out << '/' << td.all_times << '\t';
	*/
  if (td.average_nv_time > 1.0e+30)
  {
		out << "/\t/";
  }
	else
		out << td.average_nv_time << '\t' << td.reduce_time + td.average_nv_time;
  if(PRINT_HADAMARD)
    out << '\t' << td.hadamard;
  out << endl;
	return out;
}
double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}
inline double Mult(double a[DIM], double b[DIM])
{
	double ret = 0;
	for (int i = 0; i < DIM; i++)
	{
		ret += a[i] * b[i];
	}
	return ret;

}
inline void NumMult(double ret[DIM], double a[DIM], int n)
{
	for (int i = 0; i < DIM; i++)
	{
		ret[i] = a[i] * n;
	}
}
inline void AddEq(double ret[DIM], double a[DIM])
{
	for (int i = 0; i < DIM; i++)
	{
		ret[i] += a[i];
	}
}
inline void SubEq(double ret[DIM], double a[DIM])
{
	for (int i = 0; i < DIM; i++)
	{
		ret[i] -= a[i];
	}
}
int RoundToInt(double d)
{
	if (d > 0)
		return int(d + 0.5);
	else
		return int(d - 0.5);
}
//使用全局变量存储大型数组
double d_m[DIM][DIM] = { 0 };
double d_basis[DIM][DIM] = { 0 };
double d_save[DIM] = { 0 };
void NearVectorCoef(int d_coef[DIM], double d_v[DIM])
{
	//求解向量线性表示的整系数
	double mid_temp[DIM] = { 0 };
	double d_v_in[DIM];
	double d_v_out[DIM] = { 0 };
	memcpy(d_v_in, d_v, sizeof(d_v_in));
	for (int i = DIM - 1; i >= 0; i--)
	{
		int temp = RoundToInt(Mult(d_basis[i], d_v_in) / d_save[i]);
		d_coef[i] = temp;
		NumMult(mid_temp, d_m[i], temp);
		SubEq(d_v_in, mid_temp);
		AddEq(d_v_out, mid_temp);
	}
	for (int i = 0; i < DIM; i++)
	{
		noise[i] *= computed;
		computed++;
		noise[i] += (d_v_out[i] - d_v[i]) * noise_times;
		noise[i] /= computed;
	}
}
inline ZZ GetHide(int d_coef[DIM], mat_ZZ& m)
{
	ZZ d_prob_hide;
	for (int i = 1; i <= DIM; i++)
	{
		d_prob_hide += d_coef[i - 1] * m(i, DIM);
	}
	//cout << d_prob_hide << endl;
	return d_prob_hide;
}
inline ZZ GetHide(mat_ZZ& m, double d_v[DIM])
{
	int d_coef[DIM];
	NearVectorCoef(d_coef, d_v);
	return GetHide(d_coef, m);
}
void getRandList(int* list, int max_term, int choose)
{
	pair<int, int>* temp_list = new pair<int, int>[max_term];
	for (int i = 0; i < max_term; i++)
	{
		temp_list[i].first = rand();
		temp_list[i].second = i;
	}
	sort(temp_list, temp_list + max_term);
	for (int i = 0; i < choose; i++)
	{
		list[i] = temp_list[i].second;
	}
	delete[] temp_list;
}
void toString3(char* s, int n)
{
	s[0] = n % 1000 / 100 + '0';
	s[1] = n % 100 / 10 + '0';
	s[2] = n % 10 + '0';
}

char char_alpha[EQU][LEN] = { 0 };
char char_beta[EQU][LEN] = { 0 };
TestData td;

void lattice_reduction(int choose[DIM - 1], char* basis_name, int get_pid = 0)
{
	mpz_t alpha[DIM];
	for (int i = 0; i < DIM - 1; i++)
	{
		mpz_init_set_str(alpha[i], char_alpha[choose[i]], 10);
	}
	//set matrix
	ZZ_mat<mpz_t> A;
	mpz_t alpha_temp[DIM];
	for (int i = 0; i < DIM - 1; i++)
	{
		mpz_init(alpha_temp[i]);
		mpz_mul(alpha_temp[i], alpha[i], mpz_pow2[BITS]);
	}
	A.resize(DIM, DIM);
	for (int i = 0; i < DIM - 1; i++)
	{
		A[i][i] = mpz_pow2[BITS * 2];
		A[DIM - 1][i] = alpha_temp[i];
	}
	A[DIM - 1][DIM - 1] = mpz_pow2[0];
	vector<Strategy> strategies = load_strategies_json("default.json");
	double used_time=0;
	for(int i=3;i<=BKZ_BLOCK_SIZE;i++)
	{
		BKZParam b(i, strategies, 1 - (1e-10), BKZ_MAX_TIME, 0, BKZ_SET_MAX_TIME);
		auto start = high_resolution_clock::now();
		bkz_reduction(&A, NULL, b);
		auto end = high_resolution_clock::now();
		used_time += duration_cast<milliseconds>(end - start).count() / 1000.0;
	}
	td.reduce_time = used_time;

	FILE* fout;
	fout = fopen(basis_name, "w");

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			gmp_fprintf(fout, "%Zd\n", A[i][j]);
		}
	}
	fclose(fout);
}

void solve_hnp(int choose[DIM - 1], char* basis_name, int get_pid = 0)
{
	mat_ZZ m;
	vec_ZZ v;
	vec_ZZ ret;
	vec_ZZ Coeff;
	vec_ZZ Known;
	m.SetDims(DIM, DIM);
	v.SetLength(DIM);
	Coeff.SetLength(DIM);
	Known.SetLength(DIM);
	ifstream fin_basis(basis_name);
	for (int i = 1; i <= DIM; i++)
	{
		for (int j = 1; j <= DIM; j++)
		{
			fin_basis >> m(i, j);
		}
	}
	for (int i = 1; i <= DIM - 1; i++)
	{
		Coeff(i) = to_ZZ(char_alpha[choose[i - 1]]);
		Known(i) = to_ZZ(char_beta[choose[i - 1]]);
		v(i) = Known(i);
		v(i) *= ZZ_pow2[UNKNOWN] * ZZ_pow2[BITS];
		v(i) += ZZ_pow2[UNKNOWN - 1] * ZZ_pow2[BITS];
	}
	//LLL_RR(m);
	v(DIM) = 1;
	fin_basis.close();

	//计算Hadamard比率
	ZZ mo[DIM];
	for (int i = 1; i <= DIM; i++)
	{
		for (int j = 1; j <= DIM; j++)
		{
			mo[i - 1] += m(i, j) * m(i, j);
		}
		RR mid = MakeRR(mo[i - 1], 0);
		mid = sqrt(mid);
		mo[i - 1] = RoundToZZ(mid);
	}
	ZZ temp12 = to_ZZ("1");
	ZZ temp13 = to_ZZ("1");
	for (int i = 0; i < DIM - 1; i++)
		temp12 *= ZZ_pow2[BITS] * ZZ_pow2[BITS];
	for (int i = 0; i < DIM; i++)
	{
		temp13 *= mo[i];
	}
	RR temp22, temp23;
	temp22 = MakeRR(temp12, 0);
	temp23 = MakeRR(temp13, 0);
	RR temp24 = MakeRR(to_ZZ("1"), 0) / MakeRR(ZZ(DIM), 0);
  td.hadamard = pow(temp22 / temp23, temp24);
	//计算Gram-Schmidt正交基
	mat_RR m_RR;
	vec_RR v_RR;
	mat_RR basis;
	mat_RR mu;
	m_RR.SetDims(DIM, DIM);
	v_RR.SetLength(DIM);
	basis.SetDims(DIM, DIM);
	mu.SetDims(DIM, DIM);
	for (int i = 1; i <= DIM; i++)
	{
		for (int j = 1; j <= DIM; j++)
		{
			basis(i, j) = MakeRR(m(i, j), 0);
			m_RR(i, j) = MakeRR(m(i, j), 0);
		}
		v_RR(i) = MakeRR(v(i), 0);
	}
	for (int i = 2; i <= DIM; i++)
	{
		for (int j = 1; j <= i; j++)
		{
			mu(i, j) = (m_RR(i) * basis(j)) / (basis(j) * basis(j));
		}
		for (int j = 1; j <= i - 1; j++)
		{
			basis(i) -= mu(i, j) * basis(j);
		}
	}

	//求解向量线性表示方式的整系�?
	RR save[DIM + 1];
	ZZ save_coef[DIM + 1];
	for (int i = 1; i <= DIM; i++)
	{
		save[i] = basis(i) * basis(i);
	}
	//使用double类型求解可能的隐藏数
	//将基和目标向量转换为double类型
	double d_v[DIM] = { 0 };
	for (int i = 1; i <= DIM; i++)
	{
		for (int j = 1; j <= DIM; j++)
		{
			RR temp1 = MakeRR(m(i, j), 0) / RR_pow2[BITS];
			RR temp2 = basis(i, j) / RR_pow2[BITS];
			conv(d_m[i - 1][j - 1], temp1);
			conv(d_basis[i - 1][j - 1], temp2);
		}
		ZZ temp3 = v(i) / ZZ_pow2[BITS];
		conv(d_v[i - 1], temp3);
	}
	//将预计算结果转化为double类型
	for (int i = 1; i <= DIM; i++)
	{
		RR temp = save[i] / RR_pow2[2 * BITS];
		conv(d_save[i - 1], temp);
	}

	//读取正态分布数据
	const int NPD = 1024;   //正态分布数据量
	double norm_ppf[NPD] = { 0 };
	ifstream fin_np("norm_ppf_1024.txt");
	for (int i = 0; i < NPD; i++)
	{
		fin_np >> norm_ppf[i];
		norm_ppf[i] *= d_pow2[UNKNOWN - 1] / 4;
	}
	fin_np.close();

	auto start = high_resolution_clock::now();
	//ofstream fout2("r.txt");
	int count_get_hide = 0;
	//int last_get_hide=0;
	for (int z = 0; z < TEST_TIMES; z++)
	{
		double d_v_temp[DIM] = { 0 };
		memcpy(d_v_temp, d_v, sizeof(d_v_temp));
		for (int i = 0; i < DIM - 1; i++)
		{
			d_v_temp[i] += norm_ppf[rand() % NPD];
			d_v_temp[i]+=noise[i];
		}
		//d_v_temp[0] -= d_pow2_445 / 10;
		ZZ hide = GetHide(m, d_v_temp);
		//检测隐藏数是否正确
		bool get_hide = true;
		for (int i = 1; i <= DIM - 1; i++)
		{
			if ((hide * Coeff(i)) % ZZ_pow2[BITS] / ZZ_pow2[UNKNOWN] != Known(i))
			{
				get_hide = false;
				break;
			}
		}
		if (get_hide)
		{
			//cout << z << '\t' << count_get_hide<<'\t'<<hide % ZZ_pow2[BITS] << endl;
			count_get_hide++;
			//last_get_hide=z+1;
			if (count_get_hide == 1)
			{
				td.success_times = count_get_hide;
				td.all_times = z + 1;
				//cout << count_get_hide<<'/'<< z+1 << '\t';
				break;
			}
			if (PRINT_HIDE)
				cout << endl << hide % ZZ_pow2[BITS] << endl;
		}
		//fout2 << hide << endl;
		//if (z % 50000 == 0)
		//	cout << z << endl;
		if (z == TEST_TIMES - 1)
		{
			if (count_get_hide != 0)
			{
				td.success_times = count_get_hide;
				td.all_times = TEST_TIMES;
			}
			else
			{
				td.success_times = 0;
				td.all_times = TEST_TIMES;
			}
		}
	}
	//fout2.close();
	auto end = high_resolution_clock::now();
	if (count_get_hide > 0)
		td.average_nv_time = duration_cast<milliseconds>(end - start).count() / 1000.0;
	else
		td.average_nv_time = 1.0e+32;
}
int main()
{
	srand((unsigned)time(NULL));

	int process_i;
	for (process_i = 0; process_i < PROCEEDS; process_i++) {
		int pid = fork();
		if (pid == -1) {
			perror("fork error");
			exit(1);
		}
		else if (pid == 0) {
			break;
		}
	}

	for (int i = 0; i < getpid() * 11 % 2048; i++)
	{
		rand();
	}

	ifstream fin1(COEFF);
	ifstream fin2(KNOWN);
	for (int i = 0; i < EQU; i++)
	{
		fin1 >> char_alpha[i];
		fin2 >> char_beta[i];
	}
	fin1.close();
	fin2.close();
	//初始�?

	gmp_init();
	NTL_init();
	char basis_name[] = { "log/log_000_basis.txt" };
	char coeff_name[] = { "log/log_000_coeff.txt" };
	char known_name[] = { "log/log_000_known.txt" };
	toString3(basis_name + 8, getpid());
	toString3(coeff_name + 8, getpid());
	toString3(known_name + 8, getpid());
	ofstream fout1(coeff_name);
	ofstream fout2(known_name);
	//随机选择方程
	//cout<<"数据格式：\n序号\t\t约化时间\t成功率\t\t最近向量平均求解时�?<<endl;
	int begin_id = getpid() % 1000 * 100;
	for (int i = begin_id; i < begin_id + PARALLEL_TIMES; i++)
	{
		int choose[DIM - 1] = { 0 };
		getRandList(choose, EQU, DIM - 1);
		td.data_id = i;
		//����м��ļ�
		for (int j = 0; j < DIM - 1; j++)
		{
			fout1 << char_alpha[choose[j]] << endl;
			fout2 << char_beta[choose[j]] << endl;
		}

		//���Լ��
		lattice_reduction(choose, basis_name);

		//���������
		solve_hnp(choose, basis_name, i);
		cout << td;
	}
	fout1.close();
	fout2.close();
}
