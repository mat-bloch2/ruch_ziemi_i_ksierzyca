#define  _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h> 
#include <cmath>
#define dt 6000
#define maxT 30*24*3600
//zmiene i struktura do listy
int size[2] = {0,0};
struct vektor
{
	int x;
	struct vektor *next;
};
// struktury przechowujonce i segreguj¹ce dane
struct obiekt
{
	float	r[2];
	float	v[2];
	float	f[2];
	float    M;
};
struct orbita
{
	float a, c, e, f1[2], f2[2];
};
struct obiekty
{
	struct obiekt ksierzyc;
	struct obiekt ziemia;
	struct orbita orbita_ksierzyca;
	int t;
};

//funkcje pomocnicze
void zerowanie_wartosci(struct obiekty *zk);
float liczenie_odleglosci_miedzy_cialami(float *r1, float *r2);
//funkcje odpowiedzalne za uk³ad ziemia ksie¿yc
void krok_verleta_startujoncy(struct obiekt p, float *rPlus);
void metoda_numeryczna_prendkosciowa_verleta(struct obiekty *zk, FILE *plik);
void metoda_numeryczna_verleta(struct obiekt *p, float *rPlus, FILE *plik);
float energia_calkowita(struct obiekty p, float odleklosc);
void liczenie_sily(struct obiekty *zk);
void metoda_numeryczna_euler(struct obiekt *P, FILE *plik);
//funkcje odpowiedzalne za sprawdzaenie elipsy
void dane_elipsy(float *wymiar_apogeum,float *wymiar_perygeum, struct orbita *orbita);
float sprawdzenie_elipsy(struct obiekty *zk,struct vektor *Punkty);
//funkcje odpowiedzalne za liste
void push_back(struct vektor *dane, int a, int *size);
int szukaj(int i, struct vektor dane);
void destroy_vecktor(struct vektor *dane, int *size);
int  main()
{
	FILE *plik = fopen("ziemia_ksierzyc_dane.txt", "w");
	fprintf(plik, "#################ruch ziemi i ksierzyca##########################\n");
	float krPlus[2], zrPlus[2];
	struct obiekty zk;
	struct vektor lista_punktow[2];
	float  wymiar_perygeum[2],wymiar_apogeum[2],przybli¿enie_elipsy;
	zerowanie_wartosci(&zk);
	fprintf(plik, "#metoda prendkosciowa verleta: \n");
	fprintf(plik, "#t \t krx \t\t\t kry \t\t kvx \t\t kvy\t\t kfx \t\t kfy \t\tZrx \t\t Zry \t\t Zvx \t\t Zvy\t\t Zfx\t\t Zfy \t\t EC \t\t odleglosc\n");
	liczenie_sily(&zk);
	zk.orbita_ksierzyca.f1[0]= zk.ziemia.r[0];
	zk.orbita_ksierzyca.f1[1]= zk.ziemia.r[1];
	while (maxT > zk.t)
	{
		if (liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r) > 405696000)
		{
			wymiar_apogeum[0] = zk.ksierzyc.r[0];
			wymiar_apogeum[1] = zk.ksierzyc.r[1];
		}

		if (liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r) > 363104000)
		{
			wymiar_perygeum[0] = zk.ksierzyc.r[0];
			wymiar_perygeum[1] = zk.ksierzyc.r[1];
		}
		metoda_numeryczna_prendkosciowa_verleta(&zk, plik);
		push_back(&lista_punktow[0], zk.ksierzyc.r[0], &size[0]);
		push_back(&lista_punktow[1], zk.ksierzyc.r[1], &size[1]);
		//printf("%e\n", liczenie_odleglosci_miedzy_cialami(druga_ogniskowa, zk.ksierzyc.r)+liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, ogniskowa));
	}
	dane_elipsy(wymiar_apogeum, wymiar_perygeum,&zk.orbita_ksierzyca);
	przybli¿enie_elipsy=sprawdzenie_elipsy(&zk, lista_punktow);
	printf("elipsa na orbicie istnieje w przybli¿eniu ok:%e o danych a=%e\tc=%e\te=%e\tf1=(%e,%e)\tf2=(%e,%e)",przybli¿enie_elipsy ,zk.orbita_ksierzyca.a, zk.orbita_ksierzyca.c, zk.orbita_ksierzyca.e, zk.orbita_ksierzyca.f1[0], zk.orbita_ksierzyca.f1[1], zk.orbita_ksierzyca.f2[0], zk.orbita_ksierzyca.f2[1]);
	destroy_vecktor(&lista_punktow[0], &size[0]);
	destroy_vecktor(&lista_punktow[1], &size[1]);
	fprintf(plik, "\n\n");
	zerowanie_wartosci(&zk);
	fprintf(plik, "#metoda eulera:\n");
	fprintf(plik, "#t \t krx \t\t\t kry \t\t kvx \t\t kvy\t\t kfx \t\t kfy \t\tZrx \t\t Zry \t\t Zvx \t\t Zvy\t\t Zfx\t\t Zfy \t\t EC \t\t odleglosc\n");
	liczenie_sily(&zk);
	while (maxT > zk.t)
	{
		zk.t += dt;
		//printf("%e\n", liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r));
		fprintf(plik, "%d", zk.t);
		metoda_numeryczna_euler(&zk.ziemia, plik);
		metoda_numeryczna_euler(&zk.ksierzyc, plik);
		fprintf(plik, "\t%e", energia_calkowita(zk, liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r)));
		fprintf(plik, "\t%e", liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r));
		fprintf(plik, "\n");
		liczenie_sily(&zk);
	}
	fprintf(plik, "\n\n");
	zerowanie_wartosci(&zk);
	fprintf(plik, "#metoda  verleta:\n");
	fprintf(plik, "#t \t krx \t\t\t kry \t\t kvx \t\t kvy\t\t kfx \t\t kfy\t\tZrx \t\t Zry \t\t Zvx \t\t Zvy\t\t Zfx\t\t Zfy \t\t EC \t\t odleglosc\n");
	liczenie_sily(&zk);
	krok_verleta_startujoncy(zk.ziemia, zrPlus);
	krok_verleta_startujoncy(zk.ksierzyc, krPlus);
	while (maxT > zk.t)
	{
		//printf("%e\n", liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r));
		liczenie_sily(&zk);
		zk.t += dt;
		fprintf(plik, "%d", zk.t);
		metoda_numeryczna_verleta(&zk.ziemia, zrPlus, plik);
		metoda_numeryczna_verleta(&zk.ksierzyc, krPlus, plik);
		fprintf(plik, "\t%e", energia_calkowita(zk, liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r)));
		fprintf(plik, "\t%e", liczenie_odleglosci_miedzy_cialami(zk.ksierzyc.r, zk.ziemia.r));
		fprintf(plik, "\n");
	}
	fclose(plik);
	getchar();
	return 0;
}
float liczenie_odleglosci_miedzy_cialami(float *r1, float *r2)
{
	float tmp;
	tmp=pow((r1[0] - r2[0]), 2) + pow((r1[1] - r2[1]), 2);
	return sqrt(tmp);
}
void liczenie_sily(struct obiekty *zk)
{
	float odleglosc = liczenie_odleglosci_miedzy_cialami(zk->ksierzyc.r, zk->ziemia.r);
	float tmp;
	tmp = -6.67408e-11 * zk->ziemia.M * zk->ksierzyc.M / pow(odleglosc,3);
	zk->ksierzyc.f[0] = tmp * (zk->ksierzyc.r[0] - zk->ziemia.r[0]);
	zk->ksierzyc.f[1] = tmp * (zk->ksierzyc.r[1] - zk->ziemia.r[1]);
	zk->ziemia.f[0] = tmp * (zk->ziemia.r[0] - zk->ksierzyc.r[0]);
	zk->ziemia.f[1] = tmp * (zk->ziemia.r[1] - zk->ksierzyc.r[1]);
}
float energia_calkowita(struct obiekty p, float odleklosc)
{
	float EC = 0,ZEk, KEk, Ep;
	ZEk = (((p.ziemia.v[0] * p.ziemia.v[0]) + (p.ziemia.v[1] * p.ziemia.v[1])) * p.ziemia.M) / 2.0;
	KEk = (((p.ksierzyc.v[0] * p.ksierzyc.v[0]) + (p.ksierzyc.v[1] * p.ksierzyc.v[1])) * p.ksierzyc.M) / 2.0;
	Ep = p.ziemia.M*-6.67408e-11/odleklosc;
	Ep = Ep * p.ksierzyc.M;
	EC = ZEk+ KEk +Ep;
	return EC;
}
void metoda_numeryczna_euler(struct obiekt *P, FILE *plik)
{
	float EC;
	P->r[0] += P->v[0] * dt;
	P->r[1] += P->v[1] * dt;
	P->v[0] += P->f[0] / P->M * dt;
	P->v[1] += P->f[1] / P->M * dt;
	fprintf(plik, "\t %e \t %e \t %e \t %e \t %e \t %e ", P->r[0], P->r[1], P->v[0], P->v[1], P->f[0], P->f[1]);
}
void metoda_numeryczna_verleta(struct obiekt *p, float *rPlus, FILE *plik)
{
	float  rMinus[2];
	rMinus[0] = p->r[0];
	rMinus[1] = p->r[1];
	p->r[0] = rPlus[0];
	p->r[1] = rPlus[1];
	rPlus[0] = -rMinus[0] + 2 * p->r[0] + p->f[0] / p->M * dt*dt;
	rPlus[1] = -rMinus[1] + 2 * p->r[1] + p->f[1] / p->M * dt*dt;
	p->v[0] = (rPlus[0] - rMinus[0]) / (2 * dt);
	p->v[1] = (rPlus[1] - rMinus[1]) / (2 * dt);
	fprintf(plik, "\t %e \t %e \t %e \t %e \t %e \t %e ", p->r[0], p->r[1], p->v[0], p->v[1], p->f[0], p->f[1]);
}
void metoda_numeryczna_prendkosciowa_verleta(struct obiekty *zk, FILE *plik)
{
	float EC = 0, Kv_pul_t[2], Zv_pul_t[2];
	zk->t += dt;
	Kv_pul_t[0] = zk->ksierzyc.v[0] + (zk->ksierzyc.f[0] / (2 * zk->ksierzyc.M)*dt);
	Kv_pul_t[1] = zk->ksierzyc.v[1] + (zk->ksierzyc.f[1] / (2 * zk->ksierzyc.M)*dt);
	Zv_pul_t[0] = zk->ziemia.v[0] + (zk->ziemia.f[0] / (2 * zk->ziemia.M)*dt);
	Zv_pul_t[1] = zk->ziemia.v[1] + (zk->ziemia.f[1] / (2 * zk->ziemia.M)*dt);
	zk->ksierzyc.r[0] += Kv_pul_t[0] * dt;
	zk->ksierzyc.r[1] += Kv_pul_t[1] * dt;
	zk->ziemia.r[0] += Zv_pul_t[0] * dt;
	zk->ziemia.r[1] += Zv_pul_t[1] * dt;
	liczenie_sily(zk);
	zk->ksierzyc.v[0] = Kv_pul_t[0] + (zk->ksierzyc.f[0] / (2 * zk->ksierzyc.M)*dt);
	zk->ksierzyc.v[1] = Kv_pul_t[1] + (zk->ksierzyc.f[1] / (2 * zk->ksierzyc.M)*dt);
	zk->ziemia.v[0] = Zv_pul_t[0] + (zk->ziemia.f[0] / (2 * zk->ziemia.M) *dt);
	zk->ziemia.v[1] = Zv_pul_t[1] + (zk->ziemia.f[1] / (2 * zk->ziemia.M) *dt);
	EC = energia_calkowita(*zk, liczenie_odleglosci_miedzy_cialami(zk->ksierzyc.r, zk->ziemia.r));
	fprintf(plik, "%d\t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e  \n", zk->t, zk->ksierzyc.r[0], zk->ksierzyc.r[1], zk->ksierzyc.v[0], zk->ksierzyc.v[1], zk->ksierzyc.f[0], zk->ksierzyc.f[1], zk->ziemia.r[0], zk->ziemia.r[1], zk->ziemia.v[0], zk->ziemia.v[1], zk->ziemia.f[0], zk->ziemia.f[1], EC,liczenie_odleglosci_miedzy_cialami(zk->ksierzyc.r, zk->ziemia.r));
}
void zerowanie_wartosci(struct obiekty *zk)
{
	zk->t = 0;
	zk->ksierzyc.r[0] = 405696000;
	zk->ksierzyc.r[1] = 0;
	zk->ksierzyc.v[0] = 0;
	zk->ksierzyc.v[1] = 968;
	zk->ziemia.M = 5.97219e24;
	zk->ksierzyc.M = 7.347673e22;
	zk->ziemia.r[0] = 0;
	zk->ziemia.r[1] = 0;
	zk->ziemia.v[0] = 0;
	zk->ziemia.v[1] = 0;
}
void krok_verleta_startujoncy(struct obiekt p, float *rPlus)
{
	rPlus[0] = p.r[0] + p.v[0] * dt + p.f[0] / (2 * p.M) *dt*dt;
	rPlus[1] = p.r[1] + p.v[1] * dt + p.f[1] / (2 * p.M) *dt*dt;
}
void dane_elipsy( float *wymiar_apogeum, float *wymiar_perygeum, struct orbita *orbita)
{
	float srodek[2],wektor[2];
	orbita->a= liczenie_odleglosci_miedzy_cialami(wymiar_apogeum, wymiar_perygeum)/2;
	//printf("%e \t %e\n%e \t %e\t%e",wymiar_apogeum[0],wymiar_apogeum[1], wymiar_perygeum[0], wymiar_perygeum[1],a);
	srodek[0] = (wymiar_apogeum[0] + wymiar_perygeum[0]) / 2;
	srodek[1] = (wymiar_apogeum[1] + wymiar_perygeum[1]) / 2;
	orbita->c = liczenie_odleglosci_miedzy_cialami(srodek, orbita->f1);
	orbita->e = orbita->c / orbita->a;
	wektor[0] = -(orbita->f1[0] - srodek[0]);
	wektor[1] = -(orbita->f1[1] - srodek[1]);
	orbita->f2[0] = srodek[0] + wektor[0];
	orbita->f2[1] = srodek[1] + wektor[1];
	//printf("\n\t%e \t %e \n", srodek[0], srodek[1]);
	
}
void push_back(struct vektor *dane, int a, int *size)
{
	if (*size == 0)
	{
		dane->x = a;
	}
	else
	{
		struct vektor *danenext;
		danenext = dane;
		for (int i = 1; i < *size; i++)
		{
			danenext = danenext->next;
		}
		struct vektor *tmp;
		tmp = (struct vektor*)malloc(sizeof(struct vektor));
		tmp->x = a;
		tmp->next = 0;
		danenext->next = tmp;
	}
	*size = *size + 1;
}
int szukaj(int i, struct vektor dane)
{
	struct vektor *tmp = 0;
	if (i == 0)
	{
		return dane.x;
	}
	tmp = dane.next;

	for (int j = 1; j <= i; j++)
	{
		if (j == i)
			return tmp->x;
		else
			tmp = tmp->next;
	}
}
void destroy_vecktor(struct vektor *dane, int *size)
{
	int i;
	int j;
	struct vektor *danenext;
	while (*size > 0)
	{
		danenext = dane;
		for (int i = 1; i < *size - 1; i++)
		{
			danenext = danenext->next;
		}
		free(danenext->next);
		danenext->next = NULL;
		*size = *size - 1;
	}
}
float sprawdzenie_elipsy(obiekty *zk,struct vektor *Punkty)
{
	int i;
	float tmp[2], stala,max=0;
	for (i = 0; i < size[0]; i++)
	{
		tmp[0] = szukaj(i, Punkty[0]);
		tmp[1] = szukaj(i, Punkty[1]);
		stala=liczenie_odleglosci_miedzy_cialami(zk->orbita_ksierzyca.f1, tmp)+ liczenie_odleglosci_miedzy_cialami(zk->orbita_ksierzyca.f2, tmp);
		if ((max+(2 * zk->orbita_ksierzyca.a)) < stala)
		{
			max = stala - (2 * zk->orbita_ksierzyca.a);
		}
			
	}
	

	return max;
}