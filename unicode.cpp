#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>

using namespace std;

//#define DEBUG 1

//globalne potrebne za vise funckija

template <typename T>
std::string to_string(T const& value)
{
	stringstream sstr;
	sstr << value;
	return sstr.str();
}

// konstante 

double zc;

const double dZC = 10000;

const double pi = 3.14159265359;
// gravitacionalna konstanta
const double G = 6.67384E-11;
// Plankova konstanta
const double h = 6.62606957E-34;
// redukovana PK
const double hred = h / (2 * pi);
// inverzni broj elektrona po protonu
const double Y = 0.5;
// barionska masa
const double mb = 1.67262178E-27;
// brzina svetlosti
const double c = 299792458;
// masa Sunca
const double Ms = 1.9891E30;
// masa elektrona
const double Me = 9.10938291E-31;

/* ------------------------ Politropa -------------------------- */

/* Politropski indeksi */

double n;

double n_1 = 1.5; /* politropski indeks za recimo Jupiter */
double n_2 = 3.0; /* politropski indeks za wd */
double n_3 = 5.0; /* politropski indeks za globularna jata */ 

/* Promenljive teta, ksi */

double teta;
double ksi;

/*
// smene

double ro_c;
double K;
double L = 1.0/n - 1.0;
double alfa_sq = (1.0 + n) * K * pow(ro_c, L) / (4.0 * pi * G);
double alfa = pow(alfa_sq, 0.5);

// gustina

double ro = K * pow(teta, n)

*/

/* ------------------------------------------------------------- */

/* Deo ispod služi da vrednost pod korenom nikada ne bude manja od 0 */


// vraca izvod alfe po ksi

double f1(double fi0, double alfa1, double ksi)
{
	double potkorena_vrednost = fi0 * fi0 - 1 / pow(zc, 2);
	if (potkorena_vrednost < 0)
		{
			potkorena_vrednost = 0;
		}

	return -pow(potkorena_vrednost, 1.5) - 2 * alfa1 / ksi;
}


// vraća izvod alfa1 po ksi
double f2(double fi0, double eta0, double alfa2, double ksi)
{
	double potkorena_vrednost = fi0 * fi0 - 1 / pow(zc, 2);
	if (potkorena_vrednost < 0)
		{
			potkorena_vrednost = 0;
		}

	return -3 * pow(potkorena_vrednost, 0.5) * fi0 * eta0 + 1 - 2 * alfa2 / ksi;

}
// vraća izvod alfa3 po ksi
double f3(double fi0, double eta2, double alfa3, double ksi)
{
	double potkorena_vrednost = fi0 * fi0 - 1 / pow(zc, 2);
	if (potkorena_vrednost < 0)
		{
			potkorena_vrednost = 0;
		}

	return -3 * pow(potkorena_vrednost, 0.5) * fi0 * eta2 + (6. / (ksi * ksi)) * eta2 - 2 * (alfa3 / ksi);

}


/* Još konstanti */ 

const double H = pow(Me, 4);
const double I = pow(c, 5);
const double J = 3 * pow(h, 3);
double omega;
double v;
const double c1 = (pi*H*I) / J;
const double c2 = (16 * pi*pow(Me, 3) * pow(c , 3) * mb) / J;
//double rn = sqrt((2 * c1) / (pi*G)) / (c2*zc);

struct Fi0_Eta0_Eta2
{
	double Fi0, Eta0, Eta2, Ksi, Alfa1, Alfa2, Alfa3, A, v;
};

Fi0_Eta0_Eta2 FiEtaEta(double _zc, double omega, double A, ofstream& tabelaeta)
{
	Fi0_Eta0_Eta2 Vrednosti;

	zc = _zc;
	double fi0 = 1;
	double alfa1 = 0;
	double eta0 = 0;
	double alfa2 = 0;
	double eta2 = 0.00001;
	double alfa3 = 0;
	double ro_pol;
	double ro_ekv;

	A = 1 / sqrt(zc);


	/* ------------------ RK4 petlja ----------------------*/


	/* Deklarisanje početnog ksi i koraka */

	double ksi = 1e-150;
  const double dksi = 1e-3;

	while (zc * fi0 >= 1)
		{
			// povecaj korak
			

			/* ovaj rk4 racuna izvod fi0 od ksi */

			double k1 = dksi * alfa1;
			double l1 = dksi * f1(fi0, alfa1, ksi);

			double k2 = dksi * (alfa1 + l1 / 2);
			double l2 = dksi * f1(fi0 + k1 / 2, alfa1 + l1 / 2, ksi + dksi / 2);

			double k3 = dksi * (alfa1 + l2 / 2);
			double l3 = dksi * f1(fi0 + k2 / 2, alfa1 + l2 / 2, ksi + dksi / 2);

			double k4 = dksi * (alfa1 + l3);
			double l4 = dksi * f1(fi0 + k3, alfa1 + l3, ksi + dksi);

			double alfa1_diff = 1. / 6. * (l1 + 2 * l2 + 2 * l3 + l4);
			double fi01_diff = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);

			alfa1 += alfa1_diff;
			fi0 += fi01_diff;

			// ovaj rk4 racuna izvod eta0 od ksi

			k1 = dksi * alfa2;
			l1 = dksi * f2(fi0, eta0, alfa2, ksi);

			k2 = dksi * (alfa2 + l1 / 2);
			l2 = dksi * f2(fi0, eta0 + k1 / 2, alfa2 + l1 / 2, ksi + dksi / 2);

			k3 = dksi * (alfa2 + l2 / 2);
			l3 = dksi * f2(fi0, eta0 + k2 / 2, alfa2 + l2 / 2, ksi + dksi / 2);

			k4 = dksi * (alfa2 + l3);
			l4 = dksi * f2(fi0, eta0 + k3, alfa2 + l3, ksi + dksi);



			double alfa2_diff = 1. / 6. * (l1 + 2 * l2 + 2 * l3 + l4);
			double eta0_diff = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
			alfa2 += alfa2_diff;
			eta0 += eta0_diff;

			// ovaj rk4 racuna izvod eta2 od ksi

			k1 = dksi * alfa3;
			l1 = dksi * f3(fi0, eta2, alfa3, ksi);

			k2 = dksi * (alfa3 + l1 / 2);
			l2 = dksi * f3(fi0, eta2 + k1 / 2, alfa3 + l1 / 2, ksi + dksi / 2);

			k3 = dksi * (alfa3 + l2 / 2);
			l3 = dksi * f3(fi0, eta2 + k2 / 2, alfa3 + l2 / 2, ksi + dksi / 2);

			k4 = dksi * (alfa3 + l3);
			l4 = dksi * f3(fi0, eta2 + k3, alfa3 + l3, ksi + dksi);


			double alfa3_diff = 1. / 6. * (l1 + 2 * l2 + 2 * l3 + l4);
			double eta2_diff = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
			alfa3 += alfa3_diff;
			eta2 += eta2_diff;

			
      /* Definisanje parametara unutar petlje zbog ispisivanja u fajl tabelaeta.txt */
            
			A = 1 / sqrt(zc);

			double v = (omega * omega) / (2 * pi * G * c2);
      double rnn = sqrt((2 * c1) / (pi * G)) / (c2 * 10);
      double M_1 = 4 * pi * c2 * zc * zc * zc * rnn * rnn * rnn * ((-ksi * ksi * alfa1) + v * (((ksi * ksi * ksi) / 3) - (ksi * ksi * alfa2)));


			// # Deo koji će možda jednom da računa raspodelu gustine

			double fi_ekv = fi0 + v * (eta0 + (5.0/6.0) * (ksi*ksi) * eta2 * ((-0.5) / (3 * eta2 + ksi * alfa3)));
			double fi_45 = fi0 + v * (eta0 + (5.0/6.0) * (ksi*ksi) * eta2 * (sqrt(0.25) / (3 * eta2 + ksi * alfa3)));
			double fi_pol = fi0 + v * (eta0 + (5.0/6.0) * (ksi*ksi) * eta2 * (1 / (3 * eta2 + ksi * alfa3)));

			double ro_pol = c2 * zc * zc * zc * (fi_pol * pow((fi_pol - 1 / (zc * zc)), 3.0/2.0));
			double ro_45 =  c2 * zc * zc * zc * (fi_45 * pow((fi_pol - 1 / (zc * zc)), 3.0/2.0));
			double ro_ekv = c2 * zc * zc * zc * (fi_ekv * pow((fi_pol - 1 / (zc * zc)), 3.0/2.0));

			double roc = zc * zc * zc * pow((1 - 1/(zc*zc)), 1.5);

		


            /* ---------------------- For petlja koja izbacuje fajl "tabelaeta.txt" ---------------------------*/

            
			
			if(A < 5e-1 && zc * fi0 >= 1)
				{
					
            		double R_zapetlju = ksi * rnn / 1000;
					//cout << "Pa ja raaaadiiiiiiiim" << endl;
					//tabelaeta << ksi << "\t" << fi0 << "\t" << eta0 << "\t" << eta2 << "\t" << alfa1 << "\t" << alfa2 << "\t" << alfa3 << "\t" << zc << "\t" << omega << "\t"  << v << "\t" << A << "\t" << c1 << "\t" << c2 << "\t" << M_1 << endl;
            cout << R_zapetlju << "\t" << ro_ekv << "\t" << ro_pol << "\t" << ro_45 << "\t" << omega << endl;
            tabelaeta << R_zapetlju << "\t" << ro_ekv << "\t" << ro_pol << "\t" << ro_45 << endl;
				}

			
			

			// dodavanjeg novog ksi i ponovni ulazak u petlju

			ksi += dksi;

			//cout << ksi << "\t" << eta2 <<endl;
			//ispis << ksi << "\t" << eta2 <<endl;
		}

	Vrednosti.Fi0 = fi0;
	Vrednosti.Eta0 = eta0;
	Vrednosti.Eta2 = eta2;
	Vrednosti.Ksi = ksi;
	Vrednosti.Alfa1 = alfa1;
	Vrednosti.Alfa2 = alfa2;
	Vrednosti.Alfa3 = alfa3;


	return Vrednosti;
}
/* Početak strukture u kojoj se nalaze: masa, radijus, omega */

struct MasaRadijus
{
	double M, R, v, R_pol, R_ekv, uzP2, uzP0, M_nonrot, M_diff, Ksii, Eta2i, Alfa3i, Fi0;
};
//
MasaRadijus IzracunajMasuIRadijus(double zc, Fi0_Eta0_Eta2 parametri, double omega)
{
	MasaRadijus Zvezda;

	double rn = sqrt((2 * c1) / (pi * G)) / (c2 * zc);
	double v = (omega * omega) / (2 * pi * G * c2);
    double R = rn * parametri.Ksi / 1000.0;

	double M = 4 * pi * c2 * zc * zc * zc * rn * rn * rn * ((-parametri.Ksi * parametri.Ksi * parametri.Alfa1) + v * (((parametri.Ksi * parametri.Ksi * parametri.Ksi) / 3) - (parametri.Ksi * parametri.Ksi * parametri.Alfa2)));

	Zvezda.uzP0 =  parametri.Eta0;
	Zvezda.uzP2 =  - (5.0 / 6.0) * (parametri.Ksi * parametri.Ksi * parametri.Eta2) / (3 * parametri.Eta2 + parametri.Ksi  * parametri.Alfa3);
	Zvezda.R_pol = (rn / 1000.0) * (parametri.Ksi + (v / abs(parametri.Alfa1)) * (Zvezda.uzP0 + Zvezda.uzP2));
	Zvezda.R_ekv = (rn / 1000.0) * (parametri.Ksi + (v / abs(parametri.Alfa1)) * (Zvezda.uzP0 + Zvezda.uzP2 * (-0.5)));
	
	//mnozimo sa 1000 da bismo uvecali vaznost

	M = M / Ms;

	Zvezda.Fi0 = parametri.Fi0;
	
	// provera da li p2 radi

	//cout << Zvezda.Ksii << endl;
	//cout << Zvezda.uzP2 << endl;


	//Zvezda.Eta2i = parametri.Eta2;
	//Zvezda.Alfa3i = parametri.Alfa3;

	Zvezda.M = M;
	Zvezda.M_nonrot = 4*pi*rn*rn*rn*c2*zc*zc*zc*(-1)*(parametri.Ksi)*(parametri.Ksi)*(parametri.Alfa1);
 	Zvezda.M_diff = Zvezda.M - (Zvezda.M_nonrot*(1+v*(((parametri.Ksi/3)-parametri.Alfa2)/(abs(parametri.Alfa1)))))/Ms; //provera da li postoji razlika izmedju koriscenog izraza za Masu i onong kod Ch
 	Zvezda.R = R;
 	Zvezda.v = v;

 	
	return Zvezda;


/* ----------------------  GLAVNI DEO KODA -------------------------*/
}

int main()
{
    
	Fi0_Eta0_Eta2 Proba;
	MasaRadijus Zvezda;


	// double A = 1e-20;
	// double ZCinit = 1 / sqrt(A);


	ofstream autput;
	string txt = ".txt";
	autput.open("izlaz.txt");
	autput << "#" << " " << "Masa" << "\t" << "R"<< "\t" << "R_pol" << "\t" << "R_ekvator" << "\t" << "Greška_Masa" << "\t" << "Omega" << "\t" << "Zc" << "\t" <<"Fi0" << endl;
	ofstream tabelaeta;
	tabelaeta.open("tabelaeta.txt");
    //tabelaeta << "#" << " " << "ksi" << "\t" << "fi0" << "\t" << "eta0" << "\t" << "eta2" << "\t" << "alfa1" << "\t" << "alfa2" << "\t" << "alfa3" << "\t" << "zc" << "\t" << "omega" << "\t" << "v" << "\t" << "A" << "\t" << "c1" << "\t" << "c2" << endl;
    
    /* Definisanje range-a od omega */

    double omega_first = 1E-5;
    double omega_last = 1e-1;
    double d_omega = 1E-4;
    /* Definisanje range-a od A */

    double A_first = 1E-100;
    double A_last = 5e-1;
    const double dA = 0.001;

    double A = 0.01;
    omega = 0.1;

    

    //Zvezda.M = 0.56;



    /* ----------------------- Petlja koja vrti centralnu gustinu preko parametra A, i omega -------------------*/

   
	for (double omega = 0 ; omega <= omega_last; omega += d_omega)
		{

			for (double A = A_first; A < A_last; A += dA)
				{

					zc = 1 / sqrt(A);
					Proba = FiEtaEta(zc, omega, A, tabelaeta);
					Zvezda = IzracunajMasuIRadijus(zc, Proba, omega);
					double roc = zc * zc * zc * pow((1 - 1/(zc*zc)), 1.5);
					double v_rad = Zvezda.R_ekv * omega;
					double T = 2 * Zvezda.R_ekv * pi / v_rad;



					//cout  << Zvezda.M << "\t" << Zvezda.R << "\t" << Zvezda.R_pol << "\t" << Zvezda.R_ekv << "\t" << Zvezda.M_diff <<  "\t" << omega  << "\t" << zc << "\t" << Zvezda.Fi0 << "\t" << roc << endl;
					//autput << Zvezda.M << "\t" << Zvezda.R << "\t" << Zvezda.R_pol << "\t" << Zvezda.R_ekv << "\t" << Zvezda.M_diff <<  "\t" << omega  << "\t" << zc << "\t" << Zvezda.Fi0 << "\t" << roc << endl;
					
					//cout << Zvezda.M << "\t" << Zvezda.R_ekv << "\t" << omega << "\t" << v_rad << "\t" << T << endl;
					//autput << Zvezda.M << "\t" << Zvezda.R_ekv << "\t" << omega << "\t" << v_rad << "\t" << T << endl;

   						
				
					

					}
					

					
                    //double gotovo = omega/omega_last * 100;
                    //cout << gotovo << "\t" << "%" << endl;
                    
				}   

	
		}
	autput.close();
	tabelaeta.close();


	return 0;
}
