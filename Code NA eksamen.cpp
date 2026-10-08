#include <iomanip> //Indeholder funktioner til formateret input/output
#include <iostream> //Indeholder funktioner til standard input/output
#include <vector> //Indeholder funktioner til vektorer
#include <cmath> //Indeholder matematiske funktioner
#include <fstream> //Indeholder funktioner til filinput/output
#include <string> //Indeholder funktioner til strengoperationer

using namespace std;
using Vektor = vector<double>;
using Matrix = vector<Vektor>; //vector<vector<double>>

//Vigtige funktioner
int vknude(const Vektor& t, double x);
Vektor eval_spline(const Vektor& t, int k, int v, double x);
Vektor dan_knudepunkter(double tmin, double tmax, int l, int k);
Vektor bsub(const Matrix& T);
Vektor house(const Matrix& T, int j);
Matrix spline_totalmatrix(const Vektor& t, int k, const Vektor& x, const Vektor& y);
void anvend_house(Matrix& T, const Vektor& v, int j);
void houseQR(Matrix& T);

//Øvrige funktioner
void introduktion();
Vektor hent_vektor(string fn);
void udskriv_vektor(const Vektor& v);
void udskriv_matrix(const Matrix& A);
void data(double& l, double& k);
void Matrix_c (Matrix&A);

int main() {
	introduktion();

	// Indlæsning af x og y fra filer
	Vektor x = hent_vektor("x.txt");
	cout << "\nVektor x = ";
	udskriv_vektor(x);
	Vektor y = hent_vektor("y.txt");
	cout << "\nVektor y = ";
	udskriv_vektor(y);

	// Indtastning af antal intervaller og orden
	double l, k;
	data(l, k);

	//Dannelse af knudepunkter
	Vektor t = dan_knudepunkter(x[0], x[x.size()-1], l, k);
	cout << "\nVektor t = ";
	udskriv_vektor(t);

	// Dannelse af totalmatrix med De Boors algoritme
	Matrix S = spline_totalmatrix(t, k, x, y);
	Matrix_c(S);
	cout << endl << endl;
	udskriv_matrix(S);
	cout << endl << endl;

	// Løsning af kvadratminimeringsproblemet med Householders algoritme
	houseQR(S);
	Matrix_c(S);

	cout << "\n\nTotalmatricen med De Boors algoritme og Householders algortime=" << endl;
	udskriv_matrix(S);

	//Udfører baglæns substitution på total matricen og løs systemet
	Vektor c = bsub(S);

	cout << "\nKoefficienter for spline:" << endl;
	udskriv_vektor(c);

	return 0;

}

Vektor eval_spline(const Vektor& t, int k, int v, double x) { //De Boors algoritme
	Vektor dl(k-1);
	Vektor dr(k-1);
	Vektor b(k);
	double saved;

	for (int s=0; s<k-1; s++) {
		dl[s] = x-t[v-s];
		dr[s] = t[v+1+s]-x;
	}

	b[0]=1;

	// Evaluering af B-splines i x
	for (int j=1; j<k; j++) {
		saved = 0;
		for (int r=0; r<j; r++) {
			double term = b[r] / (dl[j-1-r] + dr[r]);
			b[r] = saved + (dr[r] * term);
			saved = dl[j-1-r] * term;
		}
		b[j] = saved;
	}

	return b; //Funktionsværdierne i x for de k antal B-splines, der ikke er 0, ligger nu i vektor

}

int vknude(const Vektor& t, double x){
	int v;
	const int st = t.size();
	for (v = 0; v < st-1 && !(t[v] <= x && x < t[v+1]); v++);
	return v;
}

Matrix spline_totalmatrix(const Vektor& t, int k, const Vektor& x, const Vektor& y) {
	size_t m = x.size();
	size_t n = t.size()-k;
	Matrix T(m,Vektor(n+1));

	for (size_t i=0; i<m; i++) {
		int v=vknude(t, x[i]);

		if (v>=n) v=n-1;

		// Udfylder matrix T ved at evaluere B-splines for hver x-værdi
		Vektor b = eval_spline(t, k, v, x[i]);

		for (size_t j=0; j<k; j++) {
			T[i][j+v-(k-1)]=b[j];
		}
		T[i][n] = y[i]; // Tilføj y-værdien til sidste kolonne i T
	}

	return T;
}

Vektor dan_knudepunkter(double tmin, double tmax, int l, int k) {
	int n=l+2*k-1;
	Vektor t(n);

	double step = (tmax-tmin)/l;

	//Laver de første k, laver k-1 først //ÆNDRE KOMMENTER
	for(int i=0; i<=k-1; i++){
		t[i]=tmin;
	}

	//Laver dem mellem l //ÆNDRE KOMMENTER
	for(int i=1;i<l+1;i++){
		t[i+k-1]=tmin+step*(i);
	}

	//Laver de sidste k //ÆNDRE KOMMENTER
	for(int i=l+k-1;i<n;i++){
		t[i]=tmax;
	}
	return t;
}

Vektor bsub(const Matrix& T) {
	size_t n = T[0].size()-1;
	double sum;
	Vektor c(n);

	c[n-1] = T[n-1][n]/T[n-1][n-1];
	for (int i=n-2; i>=0; i--){
		sum = 0;
		for (size_t j=i+1; j<n; j++){
			sum = sum + T[i][j]*c[j];
		}
		c[i]=((T[i][n])-sum)/(T[i][i]);
	}
	return c;

}

Vektor house(const Matrix& T, int j) {
	int m = T.size();
	double sigma = 0;
	Vektor v(m);

	for (int i=1; i<m-j; i++) {
		sigma += pow(T[j+i][j], 2);
		v[i] = T[j+i][j];
	}

	double mu = sqrt(pow(T[j][j], 2) + sigma);

	if (T[j][j] <= 0) {
		v[0] = T[j][j] - mu;
	} else {
		v[0] = -sigma / (T[j][j] + mu);
	}

	return v;
}

void anvend_house(Matrix& T, const Vektor& v, int j) {
	int m = T.size();
	int n = T[0].size()-1;
	double pi, gamma, beta;

	for (int k=0; k<=n-j; k++) {
		pi = 0;
		gamma = 0;
		for (int i=0; i<m-j; i++) {
			pi += v[i] * T[j+i][j+k];
			gamma += pow(v[i],2);
		}
		if(gamma != 0){
			beta = 2/gamma;
			for (int i=0; i<m-j; i++) {
				T[j+i][j+k] -= beta * pi * v[i];
			}
		}

	}
}

void houseQR(Matrix& T) {
	int n = T.size();

	for (int j=0; j<n-1; j++) {
		Vektor v = house(T, j);

		anvend_house(T, v, j);
	}
}

//Øvrige funktioner
void introduktion() {
	cout << "Dette program tilpasser regressionssplines." << endl;
	cout << "Til dette gøres brug af De Boors algoritme til dannelsen af totalmatrix." << endl;
	cout << "Der benyttes ligeledes Housholders algoritme til mindste kvadraters metode.\n" << endl;
}

Vektor hent_vektor(string fn) {
	Vektor v;
	ifstream fil(fn);
	if (!fil.is_open()) {
		cout << "Kunne ikke åbne filen " << fn << ".\n";
		return v;
	}
	for (double e; fil >> e;) {
		v.push_back(e);
	}
	return v;
}

void udskriv_vektor(const Vektor& v) {
	cout << "{";
	if (!v.empty()) {
		cout << v[0];
		for(size_t k = 1; k < v.size(); k++) {
			cout << ", "<< v[k];
		}
	}
	cout << "}";
}

void udskriv_matrix(const Matrix& A) {
	for (const auto& v: A) {
		for (const auto& e: v) {
			cout << setw(15) << e;
		}
		cout << endl;
	}
}

void data(double& l, double& k){
	cout << "\n\nIndtast antal intervaller: " << endl;
	cin >> l;
	cout << "Indtast orden: " << endl;
	cin >> k;
}

void Matrix_c (Matrix& A){
	for(Vektor &v:A)
		for(double &e:v)
			if(abs(e)<pow(10,-12))
				e=0;
}

Matrix mgs(Matrix& A){
	size_t n=A[0].size();
	Matrix R(n,Vektor(n));
	for (int k=0;k<n;k++){
		double nor=norm(søjle(A,k));
		R[k][k]=nor;
		for (int i=0;i<A.size();i++){
			A[i][k]=A[i][k]/nor;
		}
		for (int j=k+1;j<n;j++){
			double p=prik(søjle(A,j),søjle(A,k));
			R[k][j]=p;
			for (int l=0;l<A.size();l++){
				A[l][j]=A[l][j]-p*A[l][k];
			}
		}
	}
	return R;
}




