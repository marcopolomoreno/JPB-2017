#include <stdio.h>
#include <math.h>

#define partes 1000
#define kMax 40000
#define Pi 3.141592654

double soma,epsilon;
double gama22,gama33,gama44,gama12,gama23,gama13,gama14,gama24,gama43;
double OmegaAd[partes+1],OmegaBd[partes+1],OmegaAd0[partes+1],OmegaBd0[partes+1];
double OmegaAf[partes+1],OmegaBf[partes+1],OmegaAf0[partes+1],OmegaBf0[partes+1];
double OmegaAi[partes+1],OmegaBi[partes+1],OmegaAi0[partes+1],OmegaBi0[partes+1];
double OmegaAa[partes+1],OmegaBa[partes+1],OmegaAa0[partes+1],OmegaBa0[partes+1];
double Omegad, Omegaf, Omegaa, Omegai, Af, Ai, Ad, Aa, phif, phii, phid, phia, phiCampo,phiAzul;
double deltad, deltaf, deltaa, deltai, delta, Bf, Bi, Bd, Ba, phi13,phi14,phi12,phi23;
double AAd, BBd, AAf, BBf, AAi, BBi, AAa, BBa;
double h,t,alpha12,alpha23,alpha43,alpha14,DL,eta,s12,s23,s34,s14,mi12;
double Tp,a10,a20,a30,a40,xx,Campod,Campof,Campoi,Campoa;
int Pulsos,i,j,k,m,n,g,N,q,w,passoW,inicio,fim;
double a[17][partes+1],b[17][partes+1],k1[17][partes+1],k2[17][partes+1],k3[17][partes+1],k4[17][partes+1];

double f1(double a11, double a22, double a33, double a44,
          double a12, double b12, double a23, double b23,
          double a14, double b14, double a43, double b43,
          double a13, double b13, double a24, double b24, int i)  //sistema de 4 níveis
{
    if (i==1)  return  2*Ad*b12 - 2*Bd*a12 + 2*Aa*b14 - 2*Ba*a14 + gama22*a22 +     gama44*a44;
    if (i==2)  return -2*Ad*b12 + 2*Bd*a12 + 2*Af*b23 - 2*Bf*a23 - gama22*a22 + 0.5*gama33*a33;
    if (i==3)  return -2*Ai*b43 + 2*Bi*a43 - 2*Af*b23 + 2*Bf*a23 - gama33*a33;
    if (i==4)  return  2*Ai*b43 - 2*Bi*a43 - 2*Aa*b14 + 2*Ba*a14 - gama44*a44 + 0.5*gama33*a33;
    
    if (i==5)  return -gama12*a12 - 0.0*(deltad-deltaf)*b12 + Aa*b24 + Af*b13 - Ba*a24 - Bf*a13 - Bd*(a22-a11); //a12
    if (i==6)  return -gama12*b12 + 0.0*(deltad-deltaf)*a12 + Aa*a24 - Af*a13 + Ba*b24 - Bf*b13 + Ad*(a22-a11); //b12
    if (i==7)  return -gama23*a23 - 0.0*(deltaf-deltad)*b23 - Ad*b13 + Ai*b24 + Bd*a13 + Bi*a24 - Bf*(a33-a22); //a23
    if (i==8)  return -gama23*b23 + 0.0*(deltaf-deltad)*a23 + Ad*a13 - Ai*a24 + Bd*b13 + Bi*b24 + Af*(a33-a22); //b23
    if (i==9)  return -gama14*a14 - deltaa*b14 - Ad*b24 + Ai*b13 - Bd*a24 - Bi*a13 - Ba*(a44-a11); //a14
    if (i==10) return -gama14*b14 + deltaa*a14 + Ad*a24 - Ai*a13 - Bd*b24 - Bi*b13 + Aa*(a44-a11); //b14
    if (i==11) return -gama43*a43 - deltai*b43 - Af*b24 - Aa*b13 + Bf*a24 + Ba*a13 - Bi*(a33-a44); //a43
    if (i==12) return -gama43*b43 + deltai*a43 - Af*a24 + Aa*a13 - Bf*b24 + Ba*b13 + Ai*(a33-a44); //b43
    
    if (i==13) return -gama13*a13 - (0)*b13 - Aa*b43 + Af*b12 + Ai*b14 - Ad*b23 - Ba*a43 + Bf*a12 + Bi*a14 - Bd*a23; //a13
    if (i==14) return -gama13*b13 + (0)*a13 + Aa*a43 - Af*a12 - Ai*a14 + Ad*a23 - Ba*b43 + Bf*b12 + Bi*b14 - Bd*b23; //b13
    if (i==15) return -gama24*a24 - (deltaa-deltad)*b24 - Ad*b14 + Af*b43 + Ai*b23 - Aa*b12 + Ba*a12 - Bf*a43 - Bi*a23 + Bd*a14; //a24
    if (i==16) return -gama24*b24 + (deltaa-deltad)*a24 + Ad*a14 + Af*a43 - Ai*a23 - Aa*a12 - Ba*b12 + Bf*b43 - Bi*b23 + Bd*b14; //b24
}

main() 
{
    FILE *arquivo;
    arquivo=fopen("densidade.dat","w");
    float unidade = 1e18;
    fprintf(arquivo, "Densidade rho11 rho22 rho33 rho44 SomaRho CampoDiodo CampoAzul");fprintf(arquivo,"\n");
    fprintf(arquivo, "%.0e", unidade);fprintf(arquivo,"\n");
    
    //19-02-2014
    //Modificações: 
    //              
    //Plota Populações e coerências vs densidade
    //Considera somente um grupo de átomos, em ressonância
    //Sistema de 4 níveis com diodo, femto e campos gerados
  
    gama22=(2*Pi)*6e6;gama33=(2*Pi)*660e3;gama44=(2*Pi)*1.3e6;
    gama12=0.5*gama22;gama13=0.5*gama33;gama14=0.5*gama44;
    gama24=0.5*(gama22+gama44);gama23=0.5*(gama22+gama33);
    gama43=0.5*(gama33+gama44);mi12=3.4e5;
    
    Campod=200;s12=0.64;s23=0.81;s34=1;s14=1;
    Omegaf=1e7;Omegai=1;Omegaa=0;    //em rad/s
    passoW=1;
    inicio = 0;
    fim = 30;
    
    deltaa=0;deltai=0;deltad=0;deltaf=0;
    h=5e-12;//*10000/kMax;
    DL=1e-6*100/partes;         //passo na propagação, em m

    a10=1;                     //população inicial do estado 1
    a20=0;                     //população inicial do estado 2
    a30=0;                     //população inicial do estado 3
    a40=0;                     //população inicial do estado 4	   
	
	Omegad=Campod*sqrt(s12)*mi12;
	
    for (w=inicio;w<=fim;w++)
    {
    	t=0;	  
	  	eta=unidade*passoW*w;
	  
	  	xx=eta*0.033;  //0.033 = fator de normalização
	  	alpha12=5.3e-6*s12*xx;alpha23=4.3e-7*s23*xx;
      	alpha43=2.1e-5*s34*xx;alpha14=3.6e-7*s14*xx;   
      
      	for (q=1; q<=partes; q++)
      		a[1][q] = a10;
      
      	for (q=1; q<=partes; q++)
      		for (i=2;i<=16;i++)       
         		a[i][q] = 0;
              
        for (k=1;k<=kMax-1;k++)    //abre loop de k (temporal)
        {
        	if (k % 1000 == 0)
        		printf("*");
           	
           	for (q=1;q<=partes;q++)
           	{
          
            	if (q==1)
              	{
                	OmegaAd[q] = Omegad; OmegaBd[q] = 0;
                 	OmegaAf[q] = Omegaf; OmegaBf[q] = 0;
                 	OmegaAi[q] = Omegai; OmegaBi[q] = 0;
                 	OmegaAa[q] = Omegaa; OmegaBa[q] = 0;                 
              	}

              	Ad = OmegaAd[q]; Bd = OmegaBd[q];
              	Af = OmegaAf[q]; Bf = OmegaBf[q];
              	Ai = OmegaAi[q]; Bi = OmegaBi[q];
              	Aa = OmegaAa[q]; Ba = OmegaBa[q];             
                 
              	for (j=1;j<=16;j++)
              		k1[j][q]=f1(a[1][q],a[2][q],a[3][q],a[4][q],a[5][q],a[6][q],a[7][q],a[8][q],a[9][q],
              		a[10][q],a[11][q],a[12][q],a[13][q],a[14][q],a[15][q],a[16][q],j);
                                                           
              	for (j=1;j<=16;j++)
              		k2[j][q]=f1(a[1][q]+k1[1][q]*h/2,a[2][q]+k1[2][q]*h/2,a[3][q]+k1[3][q]*h/2,
              		a[4][q]+k1[4][q]*h/2,a[5][q]+k1[5][q]*h/2,a[6][q]+k1[6][q]*h/2,a[7][q]+k1[6][q]*h/2,
              		a[8][q]+k1[8][q]*h/2,a[9][q]+k1[9][q]*h/2,a[10][q]+k1[10][q]*h/2,a[11][q]+k1[11][q]*h/2,
              		a[12][q]+k1[12][q]*h/2,a[13][q]+k1[13][q]*h/2,a[14][q]+k1[14][q]*h/2,a[15][q]+k1[15][q]*h/2,
              		a[16][q]+k1[16][q]*h/2,j);
                  
              	for (j=1;j<=16;j++)
              		k3[j][q]=f1(a[1][q]+k2[1][q]*h/2,a[2][q]+k2[2][q]*h/2,a[3][q]+k2[3][q]*h/2,
              		a[4][q]+k2[4][q]*h/2,a[5][q]+k2[5][q]*h/2,a[6][q]+k2[6][q]*h/2,a[7][q]+k2[7][q]*h/2,
              		a[8][q]+k2[8][q]*h/2,a[9][q]+k2[9][q]*h/2,a[10][q]+k2[10][q]*h/2,a[11][q]+k2[11][q]*h/2,
              		a[12][q]+k2[12][q]*h/2,a[13][q]+k2[13][q]*h/2,a[14][q]+k2[14][q]*h/2,a[15][q]+k2[15][q]*h/2,
              		a[16][q]+k2[16][q]*h/2,j);
                  
              	for (j=1;j<=16;j++)           
              		k4[j][q]=f1(a[1][q]+k3[1][q]*h,a[2][q]+k3[2][q]*h,a[3][q]+k3[3][q]*h,a[4][q]+k3[4][q]*h,
              		a[5][q]+k3[5][q]*h,a[6][q]+k3[6][q]*h,a[7][q]+k3[7][q]*h,a[8][q]+k3[8][q]*h,
              		a[9][q]+k3[9][q]*h,a[10][q]+k3[10][q]*h,a[11][q]+k3[11][q]*h,a[12][q]+k3[12][q]*h,a[13][q]+k3[13][q]*h,
              		a[14][q]+k3[14][q]*h,a[15][q]+k3[15][q]*h,a[16][q]+k3[16][q]*h,j);

              	for (j=1;j<=16;j++)
              		b[j][q]=a[j][q]+h*(k1[j][q]/6+k2[j][q]/3+k3[j][q]/3+k4[j][q]/6);   
                 
              	for (m=1;m<=16;m++)
              		a[m][q]=b[m][q];   
            
              	//equações da propagação
              	AAd = Ad + b[6][q] *alpha12*DL; BBd = Bd - b[5][q]* alpha12*DL;
              	AAf = Af + b[8][q] *alpha23*DL; BBf = Bf - b[7][q]* alpha23*DL;
              	AAa = Aa + b[10][q]*alpha14*DL; BBa = Ba - b[9][q]* alpha14*DL;
              	AAi = Ai + b[12][q]*alpha43*DL; BBi = Bi - b[11][q]*alpha43*DL;
			  
			  	if (q==partes && k==kMax-1)
			  	{
			  		soma=b[1][q]+b[2][q]+b[3][q]+b[4][q];
              		printf("\n");
      		  		printf("%d %10.8f %10.8f %10.8f %10.8f %12.10f %12.10f", 
                		w,b[1][q],b[2][q],b[3][q],b[4][q],soma,sqrt(AAa*AAa+BBa*BBa));
              		printf("\n");
              		fprintf(arquivo,"%d %12.10f %12.10f %12.10f %12.10f %14.12f %16.14f %16.14f",
                		w*passoW,b[1][q],b[2][q],b[3][q],b[4][q],soma,sqrt(AAd*AAd+BBd*BBd)/Omegad,sqrt(AAa*AAa+BBa*BBa));   
              		fprintf(arquivo,"\n");
			  		t=t+h;
		      	}
                                    
              	OmegaAd[q+1] = AAd; OmegaBd[q+1] = BBd;
              	OmegaAf[q+1] = AAf; OmegaBf[q+1] = BBf;
              	OmegaAi[q+1] = AAi; OmegaBi[q+1] = BBi;
              	OmegaAa[q+1] = AAa; OmegaBa[q+1] = BBa;            
            
           	}//fecha loop de q
       
        }//fecha loop de k

    }
    
   fclose(arquivo);
   printf("\a");
}
