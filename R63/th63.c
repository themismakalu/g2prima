#include <stdio.h>
#include <time.h>
int		ho,kcc,mi,se,npar,op,orp,r=63,rg,vp,cc[1003],dd[1003],f[1003],g[1003],nf[14];
unsigned long 	ce,c,dp,e1,i,j,k,lex,pares[546000][4],qmag[1003][1003];
double 		colump,fila;
char		nomfich[14]={'t','h'},letras[10]={'0','1','2','3','4','5','6','7','8','9'};

struct tm* tptr;
FILE* F1;		/*	themis :  compone matriz de primos concentrica de rg n  por tramos de rango     */
main()
{
    F1=fopen("tapas","r");

    fscanf(F1,"%d %lu",&rg,&ce);
    qmag[(rg-1)/2][(rg-1)/2]=ce;

    for(k=0;k>=0;k++)
	{   fscanf(F1,"%lu",&pares[k][0]);
		if(pares[k][0]==2)  break;
		fscanf(F1,"%lu",&pares[k][1]);   }
    npar=k;
    printf("\n  npar:%d\n",npar);
    fclose(F1);
    
    reortapas();

     while(r>1) 
  {	
	hora();	ho=tptr->tm_hour;	mi=tptr->tm_min;	se=tptr->tm_sec;
  	ordenpares();

	if(op==0)
  	{   printf("   matriz de rango %d completada",r);
	    hora();    printf(" %2dh %2d\' %2d\"\n",tptr->tm_hour, tptr->tm_min, tptr->tm_sec); 

	    if(r==3)   {   entacade();
    			   F1=fopen(nomfich,"w");
			   escribir();
			   printf("\n  Generada matriz de rango %d\n",rg);	}	}

   r-=2;	 }     	 /*      fin  while  r      */ 

   printf("\n finalizado a las");	hora();	   printf(" %2dh %2d\' %2d\"\n",tptr->tm_hour, tptr->tm_min, tptr->tm_sec); 
    }      /*   fin programa    */ 



reortapas()
{  F1=fopen("tapas","w");
   fprintf(F1," %d  %lu\n",rg,ce);
   for(k=1;k<npar;k++)
   {	fprintf(F1," %8lu %8lu   ",pares[k][0],pares[k][1]);
	if(k/5*5==k)  fprintf(F1,"\n");		}

	fprintf(F1," %8lu %8lu\n 2",pares[0][0],pares[0][1]);
	fclose (F1);
   }


entacade()
{
	j=13;   e1=qmag[0][0];

    while(e1>9)
   {	j--;
	nf[j]=e1%10;
	e1/=10;		 }
    nf[j-1]=e1;

    for(i=2;i<16-j;i++)  nf[i]=nf[j-3+i];

    for(e1=2;e1<16-j;e1++)  nomfich[e1]=letras[nf[e1]];

    nomfich[e1+1]='\0';
  }



ajupares()
{ 
	dp=-1;
	for(c=0;c<r-1;c++)   pares[cc[c]][0]=0;
	for(c=0;c<r-2;c++)   pares[f[c]][0]=0;
	pares[lex][0]=0;
	for(c=0;c<npar;c++)
	{   if(pares[c][0]==0)  continue; 
	    dp++;
	    pares[dp][2]=pares[c][0];	pares[dp][3]=pares[c][1];	}
	
	npar-=2*r-2;
	dp=0;
	orp+=npar/5;   if(orp > npar)  orp=0;
	for(c=orp;c<npar;c++)	{  pares[dp][0]=pares[c][2];	 pares[dp][1]=pares[c][3];    dp++; }
	for(c=0;c<orp;c++)	{  pares[dp][0]=pares[c][2];	 pares[dp][1]=pares[c][3];    dp++; }
   }




filas()
{
   fila=0;
   for(k=0;k<r-1;k+=2)   fila+=pares[cc[k]][dd[k]]+2*ce-pares[cc[k+1]][dd[k+1]];
   qmag[(rg-r)/2][(rg+r-2)/2]=r*ce-fila;
   vp=1;
   for(k=0;k<npar;k++)
       {       if(qmag[(rg-r)/2][(rg+r-2)/2]==pares[k][0] || qmag[(rg-r)/2][(rg+r-2)/2]==pares[k][1])    {      lex=k;   vp=0;   break;     }   }
   if(vp==0)  
      {     for(k=0;k<=r-2;k++) 
      {     if(qmag[(rg-r)/2][(rg+r-2)/2]!=pares[cc[k]][0] && qmag[(rg-r)/2][(rg+r-2)/2]!=pares[cc[k]][1]) continue;          vp=1;  break;      }       }
   if(vp==0) 
   {	qmag[(rg+r-2)/2][(rg-r)/2]=2*ce-qmag[(rg-r)/2][(rg+r-2)/2];
	colum();	}
           }



implementa()
   {
     qmag[(rg-r)/2][(rg-r)/2]=pares[cc[0]][dd[0]];                qmag[(rg+r-2)/2][(rg+r-2)/2]=2*ce-pares[cc[0]][dd[0]]; 
     qmag[(rg-r)/2][(rg-r)/2+r-2]=2*ce-pares[cc[r-2]][dd[r-2]];   qmag[(rg+r-2)/2][(rg-r)/2+r-2]=pares[cc[r-2]][dd[r-2]]; 
     							    
     for(k=1;k<r-2;k+=2)
    {   qmag[(rg-r)/2][(rg-r)/2+k]=2*ce-pares[cc[k]][dd[k]]; 	qmag[(rg-r)/2][(rg-r)/2+k+1]=pares[cc[k+1]][dd[k+1]];
        qmag[(rg+r-2)/2][(rg-r)/2+k]=pares[cc[k]][dd[k]];	qmag[(rg+r-2)/2][(rg-r)/2+k+1]=2*ce-pares[cc[k+1]][dd[k+1]];	}

     qmag[(rg+r-2)/2-1][(rg-r)/2]=2*ce-pares[f[r-3]][g[r-3]];	qmag[(rg+r-2)/2-1][(rg+r-2)/2]=pares[f[r-3]][g[r-3]]; 

     for(k=1;k<r-2;k+=2)
    {   qmag[(rg-r)/2+k][(rg-r)/2]=2*ce-pares[f[k-1]][g[k-1]]; 	qmag[(rg-r)/2+k][(rg+r-2)/2]=pares[f[k-1]][g[k-1]];  
        qmag[(rg-r)/2+k+1][(rg-r)/2]=pares[f[k]][g[k]];		qmag[(rg-r)/2+k+1][(rg+r-2)/2]=2*ce-pares[f[k]][g[k]];	}

   }



hora()	{	time_t t = time(NULL);	   tptr = localtime(&t);      }



escribir()
{
	hora();

	for(k=0;k<rg;k++)
	{	for(j=0;j<rg;j++)   fprintf(F1," %8lu",qmag[k][j]);   fprintf(F1,"\n");	}

   fclose(F1);
    }





ordenpares()
{
    op=1;
    for(cc[0]=0;cc[0]<npar;cc[0]++)
{   for(dd[0]=0;dd[0]<=1;dd[0]++)
{    for(cc[1]=0;cc[1]<npar; cc[1]++)  {    if(cc[1]==cc[0])   continue;
    for(dd[1]=0;dd[1]<=1;dd[1]++)
{     if(r==3)  {   filas();    if(vp==1)  continue;
                     if(op==0)  {   implementa();    ajupares();    break;  }   else  continue;   }

    for(cc[2]=0;cc[2]<npar; cc[2]++)
{      kcc=0;     for(k=0;k<2;k++)  {   if(cc[2]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[2]=0;dd[2]<=1;dd[2]++)

{   for(cc[3]=0;cc[3]<npar; cc[3]++)
{      kcc=0;     for(k=0;k<3;k++)  {   if(cc[3]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[3]=0;dd[3]<=1;dd[3]++)

{    if(r==5)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[4]=0;cc[4]<npar; cc[4]++)
{      kcc=0;     for(k=0;k<4;k++)  {   if(cc[4]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[4]=0;dd[4]<=1;dd[4]++)

{   for(cc[5]=0;cc[5]<npar; cc[5]++)
{      kcc=0;     for(k=0;k<5;k++)  {   if(cc[5]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[5]=0;dd[5]<=1;dd[5]++)

{    if(r==7)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[6]=0;cc[6]<npar; cc[6]++)
{      kcc=0;     for(k=0;k<6;k++)  {   if(cc[6]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[6]=0;dd[6]<=1;dd[6]++)

{   for(cc[7]=0;cc[7]<npar; cc[7]++)
{      kcc=0;     for(k=0;k<7;k++)  {   if(cc[7]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[7]=0;dd[7]<=1;dd[7]++)

{    if(r==9)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[8]=0;cc[8]<npar; cc[8]++)
{      kcc=0;     for(k=0;k<8;k++)  {   if(cc[8]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[8]=0;dd[8]<=1;dd[8]++)

{   for(cc[9]=0;cc[9]<npar; cc[9]++)
{      kcc=0;     for(k=0;k<9;k++)  {   if(cc[9]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[9]=0;dd[9]<=1;dd[9]++)

{    if(r==11)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[10]=0;cc[10]<npar; cc[10]++)
{      kcc=0;     for(k=0;k<10;k++)  {   if(cc[10]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[10]=0;dd[10]<=1;dd[10]++)

{   for(cc[11]=0;cc[11]<npar; cc[11]++)
{      kcc=0;     for(k=0;k<11;k++)  {   if(cc[11]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[11]=0;dd[11]<=1;dd[11]++)

{    if(r==13)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[12]=0;cc[12]<npar; cc[12]++)
{      kcc=0;     for(k=0;k<12;k++)  {   if(cc[12]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[12]=0;dd[12]<=1;dd[12]++)

{   for(cc[13]=0;cc[13]<npar; cc[13]++)
{      kcc=0;     for(k=0;k<13;k++)  {   if(cc[13]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[13]=0;dd[13]<=1;dd[13]++)

{    if(r==15)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[14]=0;cc[14]<npar; cc[14]++)
{      kcc=0;     for(k=0;k<14;k++)  {   if(cc[14]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[14]=0;dd[14]<=1;dd[14]++)

{   for(cc[15]=0;cc[15]<npar; cc[15]++)
{      kcc=0;     for(k=0;k<15;k++)  {   if(cc[15]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[15]=0;dd[15]<=1;dd[15]++)

{    if(r==17)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[16]=0;cc[16]<npar; cc[16]++)
{      kcc=0;     for(k=0;k<16;k++)  {   if(cc[16]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[16]=0;dd[16]<=1;dd[16]++)

{   for(cc[17]=0;cc[17]<npar; cc[17]++)
{      kcc=0;     for(k=0;k<17;k++)  {   if(cc[17]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[17]=0;dd[17]<=1;dd[17]++)

{    if(r==19)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[18]=0;cc[18]<npar; cc[18]++)
{      kcc=0;     for(k=0;k<18;k++)  {   if(cc[18]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[18]=0;dd[18]<=1;dd[18]++)

{   for(cc[19]=0;cc[19]<npar; cc[19]++)
{      kcc=0;     for(k=0;k<19;k++)  {   if(cc[19]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[19]=0;dd[19]<=1;dd[19]++)

{    if(r==21)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[20]=0;cc[20]<npar; cc[20]++)
{      kcc=0;     for(k=0;k<20;k++)  {   if(cc[20]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[20]=0;dd[20]<=1;dd[20]++)

{   for(cc[21]=0;cc[21]<npar; cc[21]++)
{      kcc=0;     for(k=0;k<21;k++)  {   if(cc[21]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[21]=0;dd[21]<=1;dd[21]++)

{    if(r==23)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[22]=0;cc[22]<npar; cc[22]++)
{      kcc=0;     for(k=0;k<22;k++)  {   if(cc[22]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[22]=0;dd[22]<=1;dd[22]++)

{   for(cc[23]=0;cc[23]<npar; cc[23]++)
{      kcc=0;     for(k=0;k<23;k++)  {   if(cc[23]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[23]=0;dd[23]<=1;dd[23]++)

{    if(r==25)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[24]=0;cc[24]<npar; cc[24]++)
{      kcc=0;     for(k=0;k<24;k++)  {   if(cc[24]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[24]=0;dd[24]<=1;dd[24]++)

{   for(cc[25]=0;cc[25]<npar; cc[25]++)
{      kcc=0;     for(k=0;k<25;k++)  {   if(cc[25]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[25]=0;dd[25]<=1;dd[25]++)

{    if(r==27)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[26]=0;cc[26]<npar; cc[26]++)
{      kcc=0;     for(k=0;k<26;k++)  {   if(cc[26]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[26]=0;dd[26]<=1;dd[26]++)

{   for(cc[27]=0;cc[27]<npar; cc[27]++)
{      kcc=0;     for(k=0;k<27;k++)  {   if(cc[27]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[27]=0;dd[27]<=1;dd[27]++)

{    if(r==29)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[28]=0;cc[28]<npar; cc[28]++)
{      kcc=0;     for(k=0;k<28;k++)  {   if(cc[28]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[28]=0;dd[28]<=1;dd[28]++)

{   for(cc[29]=0;cc[29]<npar; cc[29]++)
{      kcc=0;     for(k=0;k<29;k++)  {   if(cc[29]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[29]=0;dd[29]<=1;dd[29]++)

{    if(r==31)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[30]=0;cc[30]<npar; cc[30]++)
{      kcc=0;     for(k=0;k<30;k++)  {   if(cc[30]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[30]=0;dd[30]<=1;dd[30]++)

{   for(cc[31]=0;cc[31]<npar; cc[31]++)
{      kcc=0;     for(k=0;k<31;k++)  {   if(cc[31]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[31]=0;dd[31]<=1;dd[31]++)

{    if(r==33)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[32]=0;cc[32]<npar; cc[32]++)
{      kcc=0;     for(k=0;k<32;k++)  {   if(cc[32]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[32]=0;dd[32]<=1;dd[32]++)

{   for(cc[33]=0;cc[33]<npar; cc[33]++)
{      kcc=0;     for(k=0;k<33;k++)  {   if(cc[33]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[33]=0;dd[33]<=1;dd[33]++)

{    if(r==35)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[34]=0;cc[34]<npar; cc[34]++)
{      kcc=0;     for(k=0;k<34;k++)  {   if(cc[34]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[34]=0;dd[34]<=1;dd[34]++)

{   for(cc[35]=0;cc[35]<npar; cc[35]++)
{      kcc=0;     for(k=0;k<35;k++)  {   if(cc[35]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[35]=0;dd[35]<=1;dd[35]++)

{    if(r==37)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[36]=0;cc[36]<npar; cc[36]++)
{      kcc=0;     for(k=0;k<36;k++)  {   if(cc[36]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[36]=0;dd[36]<=1;dd[36]++)

{   for(cc[37]=0;cc[37]<npar; cc[37]++)
{      kcc=0;     for(k=0;k<37;k++)  {   if(cc[37]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[37]=0;dd[37]<=1;dd[37]++)

{    if(r==39)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[38]=0;cc[38]<npar; cc[38]++)
{      kcc=0;     for(k=0;k<38;k++)  {   if(cc[38]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[38]=0;dd[38]<=1;dd[38]++)

{   for(cc[39]=0;cc[39]<npar; cc[39]++)
{      kcc=0;     for(k=0;k<39;k++)  {   if(cc[39]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[39]=0;dd[39]<=1;dd[39]++)

{    if(r==41)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[40]=0;cc[40]<npar; cc[40]++)
{      kcc=0;     for(k=0;k<40;k++)  {   if(cc[40]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[40]=0;dd[40]<=1;dd[40]++)

{   for(cc[41]=0;cc[41]<npar; cc[41]++)
{      kcc=0;     for(k=0;k<41;k++)  {   if(cc[41]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[41]=0;dd[41]<=1;dd[41]++)

{    if(r==43)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[42]=0;cc[42]<npar; cc[42]++)
{      kcc=0;     for(k=0;k<42;k++)  {   if(cc[42]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[42]=0;dd[42]<=1;dd[42]++)

{   for(cc[43]=0;cc[43]<npar; cc[43]++)
{      kcc=0;     for(k=0;k<43;k++)  {   if(cc[43]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[43]=0;dd[43]<=1;dd[43]++)

{    if(r==45)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[44]=0;cc[44]<npar; cc[44]++)
{      kcc=0;     for(k=0;k<44;k++)  {   if(cc[44]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[44]=0;dd[44]<=1;dd[44]++)

{   for(cc[45]=0;cc[45]<npar; cc[45]++)
{      kcc=0;     for(k=0;k<45;k++)  {   if(cc[45]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[45]=0;dd[45]<=1;dd[45]++)

{    if(r==47)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[46]=0;cc[46]<npar; cc[46]++)
{      kcc=0;     for(k=0;k<46;k++)  {   if(cc[46]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[46]=0;dd[46]<=1;dd[46]++)

{   for(cc[47]=0;cc[47]<npar; cc[47]++)
{      kcc=0;     for(k=0;k<47;k++)  {   if(cc[47]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[47]=0;dd[47]<=1;dd[47]++)

{    if(r==49)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[48]=0;cc[48]<npar; cc[48]++)
{      kcc=0;     for(k=0;k<48;k++)  {   if(cc[48]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[48]=0;dd[48]<=1;dd[48]++)

{   for(cc[49]=0;cc[49]<npar; cc[49]++)
{      kcc=0;     for(k=0;k<49;k++)  {   if(cc[49]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[49]=0;dd[49]<=1;dd[49]++)

{    if(r==51)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[50]=0;cc[50]<npar; cc[50]++)
{      kcc=0;     for(k=0;k<50;k++)  {   if(cc[50]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[50]=0;dd[50]<=1;dd[50]++)

{   for(cc[51]=0;cc[51]<npar; cc[51]++)
{      kcc=0;     for(k=0;k<51;k++)  {   if(cc[51]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[51]=0;dd[51]<=1;dd[51]++)

{    if(r==53)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[52]=0;cc[52]<npar; cc[52]++)
{      kcc=0;     for(k=0;k<52;k++)  {   if(cc[52]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[52]=0;dd[52]<=1;dd[52]++)

{   for(cc[53]=0;cc[53]<npar; cc[53]++)
{      kcc=0;     for(k=0;k<53;k++)  {   if(cc[53]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[53]=0;dd[53]<=1;dd[53]++)

{    if(r==55)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[54]=0;cc[54]<npar; cc[54]++)
{      kcc=0;     for(k=0;k<54;k++)  {   if(cc[54]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[54]=0;dd[54]<=1;dd[54]++)

{   for(cc[55]=0;cc[55]<npar; cc[55]++)
{      kcc=0;     for(k=0;k<55;k++)  {   if(cc[55]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[55]=0;dd[55]<=1;dd[55]++)

{    if(r==57)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[56]=0;cc[56]<npar; cc[56]++)
{      kcc=0;     for(k=0;k<56;k++)  {   if(cc[56]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[56]=0;dd[56]<=1;dd[56]++)

{   for(cc[57]=0;cc[57]<npar; cc[57]++)
{      kcc=0;     for(k=0;k<57;k++)  {   if(cc[57]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[57]=0;dd[57]<=1;dd[57]++)

{    if(r==59)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[58]=0;cc[58]<npar; cc[58]++)
{      kcc=0;     for(k=0;k<58;k++)  {   if(cc[58]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[58]=0;dd[58]<=1;dd[58]++)

{   for(cc[59]=0;cc[59]<npar; cc[59]++)
{      kcc=0;     for(k=0;k<59;k++)  {   if(cc[59]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[59]=0;dd[59]<=1;dd[59]++)

{    if(r==61)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

    for(cc[60]=0;cc[60]<npar; cc[60]++)
{      kcc=0;     for(k=0;k<60;k++)  {   if(cc[60]==cc[k])  {  kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[60]=0;dd[60]<=1;dd[60]++)

{   for(cc[61]=0;cc[61]<npar; cc[61]++)
{      kcc=0;     for(k=0;k<61;k++)  {   if(cc[61]==cc[k])  {   kcc=1;   break;  }   }     if(kcc==1)  continue;
    for(dd[61]=0;dd[61]<=1;dd[61]++)

{    if(r==63)  {   filas();     if(vp==1)  continue;
                    if (op==0)  {   implementa();    ajupares();    break;  }    }

   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }   if(op==0)  break; }
   } 



  colum()
{ 
    for(f[0]=0;f[0]<npar;f[0]++)
{     if(pares[f[0]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[0]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
      kcc=0;     for(k=0;k<=r-2; k++)  {   if(f[0] == cc[k])  {   kcc=1;    break;    }   }    if(kcc==1)   continue;
    for(g[0]=0;g[0]<=1;g[0]++)

{     if(r==3)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;     continue;   }
          op=0;    break;   }

    for(f[1]=0;f[1]<npar;f[1]++)
{      kcc=0;    for(k=0;k<1;k++) {  if(f[1]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[1]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[1]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[1] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[1]=0;g[1]<=1;g[1]++)

{   for(f[2]=0;f[2]<npar;f[2]++)
{      kcc=0;    for(k=0;k<2;k++)  {  if(f[2]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[2]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[2]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[2] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[2]=0;g[2]<=1;g[2]++)

{      if(r==5)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[3]=0;f[3]<npar;f[3]++)
{      kcc=0;    for(k=0;k<3;k++) {  if(f[3]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[3]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[3]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[3] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[3]=0;g[3]<=1;g[3]++)

{   for(f[4]=0;f[4]<npar;f[4]++)
{      kcc=0;    for(k=0;k<4;k++)  {  if(f[4]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[4]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[4]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[4] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[4]=0;g[4]<=1;g[4]++)

{      if(r==7)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[5]=0;f[5]<npar;f[5]++)
{      kcc=0;    for(k=0;k<5;k++) {  if(f[5]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[5]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[5]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[5] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[5]=0;g[5]<=1;g[5]++)

{   for(f[6]=0;f[6]<npar;f[6]++)
{      kcc=0;    for(k=0;k<6;k++)  {  if(f[6]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[6]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[6]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[6] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[6]=0;g[6]<=1;g[6]++)

{      if(r==9)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[7]=0;f[7]<npar;f[7]++)
{      kcc=0;    for(k=0;k<7;k++) {  if(f[7]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[7]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[7]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[7] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[7]=0;g[7]<=1;g[7]++)

{   for(f[8]=0;f[8]<npar;f[8]++)
{      kcc=0;    for(k=0;k<8;k++)  {  if(f[8]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[8]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[8]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[8] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[8]=0;g[8]<=1;g[8]++)

{      if(r==11)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[9]=0;f[9]<npar;f[9]++)
{      kcc=0;    for(k=0;k<9;k++) {  if(f[9]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[9]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[9]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[9] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[9]=0;g[9]<=1;g[9]++)

{   for(f[10]=0;f[10]<npar;f[10]++)
{      kcc=0;    for(k=0;k<10;k++)  {  if(f[10]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[10]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[10]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[10] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[10]=0;g[10]<=1;g[10]++)

{      if(r==13)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[11]=0;f[11]<npar;f[11]++)
{      kcc=0;    for(k=0;k<11;k++) {  if(f[11]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[11]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[11]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[11] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[11]=0;g[11]<=1;g[11]++)

{   for(f[12]=0;f[12]<npar;f[12]++)
{      kcc=0;    for(k=0;k<12;k++)  {  if(f[12]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[12]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[12]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[12] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[12]=0;g[12]<=1;g[12]++)

{      if(r==15)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[13]=0;f[13]<npar;f[13]++)
{      kcc=0;    for(k=0;k<13;k++) {  if(f[13]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[13]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[13]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[13] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[13]=0;g[13]<=1;g[13]++)

{   for(f[14]=0;f[14]<npar;f[14]++)
{      kcc=0;    for(k=0;k<14;k++)  {  if(f[14]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[14]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[14]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[14] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[14]=0;g[14]<=1;g[14]++)

{      if(r==17)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[15]=0;f[15]<npar;f[15]++)
{      kcc=0;    for(k=0;k<15;k++) {  if(f[15]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[15]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[15]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[15] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[15]=0;g[15]<=1;g[15]++)

{   for(f[16]=0;f[16]<npar;f[16]++)
{      kcc=0;    for(k=0;k<16;k++)  {  if(f[16]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[16]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[16]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[16] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[16]=0;g[16]<=1;g[16]++)

{      if(r==19)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[17]=0;f[17]<npar;f[17]++)
{      kcc=0;    for(k=0;k<17;k++) {  if(f[17]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[17]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[17]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[17] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[17]=0;g[17]<=1;g[17]++)

{   for(f[18]=0;f[18]<npar;f[18]++)
{      kcc=0;    for(k=0;k<18;k++)  {  if(f[18]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[18]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[18]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[18] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[18]=0;g[18]<=1;g[18]++)

{      if(r==21)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[19]=0;f[19]<npar;f[19]++)
{      kcc=0;    for(k=0;k<19;k++) {  if(f[19]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[19]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[19]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[19] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[19]=0;g[19]<=1;g[19]++)

{   for(f[20]=0;f[20]<npar;f[20]++)
{      kcc=0;    for(k=0;k<20;k++)  {  if(f[20]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[20]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[20]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[20] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[20]=0;g[20]<=1;g[20]++)

{      if(r==23)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[21]=0;f[21]<npar;f[21]++)
{      kcc=0;    for(k=0;k<21;k++) {  if(f[21]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[21]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[21]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[21] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[21]=0;g[21]<=1;g[21]++)

{   for(f[22]=0;f[22]<npar;f[22]++)
{      kcc=0;    for(k=0;k<22;k++)  {  if(f[22]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[22]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[22]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[22] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[22]=0;g[22]<=1;g[22]++)

{      if(r==25)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[23]=0;f[23]<npar;f[23]++)
{      kcc=0;    for(k=0;k<23;k++) {  if(f[23]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[23]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[23]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[23] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[23]=0;g[23]<=1;g[23]++)

{   for(f[24]=0;f[24]<npar;f[24]++)
{      kcc=0;    for(k=0;k<24;k++)  {  if(f[24]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[24]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[24]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[24] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[24]=0;g[24]<=1;g[24]++)

{      if(r==27)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[25]=0;f[25]<npar;f[25]++)
{      kcc=0;    for(k=0;k<25;k++) {  if(f[25]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[25]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[25]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[25] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[25]=0;g[25]<=1;g[25]++)

{   for(f[26]=0;f[26]<npar;f[26]++)
{      kcc=0;    for(k=0;k<26;k++)  {  if(f[26]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[26]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[26]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[26] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[26]=0;g[26]<=1;g[26]++)

{      if(r==29)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[27]=0;f[27]<npar;f[27]++)
{      kcc=0;    for(k=0;k<27;k++) {  if(f[27]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[27]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[27]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[27] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[27]=0;g[27]<=1;g[27]++)

{   for(f[28]=0;f[28]<npar;f[28]++)
{      kcc=0;    for(k=0;k<28;k++)  {  if(f[28]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[28]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[28]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[28] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[28]=0;g[28]<=1;g[28]++)

{      if(r==31)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[29]=0;f[29]<npar;f[29]++)
{      kcc=0;    for(k=0;k<29;k++) {  if(f[29]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[29]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[29]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[29] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[29]=0;g[29]<=1;g[29]++)

{   for(f[30]=0;f[30]<npar;f[30]++)
{      kcc=0;    for(k=0;k<30;k++)  {  if(f[30]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[30]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[30]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[30] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[30]=0;g[30]<=1;g[30]++)

{      if(r==33)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[31]=0;f[31]<npar;f[31]++)
{      kcc=0;    for(k=0;k<31;k++) {  if(f[31]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[31]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[31]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[31] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[31]=0;g[31]<=1;g[31]++)

{   for(f[32]=0;f[32]<npar;f[32]++)
{      kcc=0;    for(k=0;k<32;k++)  {  if(f[32]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[32]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[32]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[32] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[32]=0;g[32]<=1;g[32]++)

{      if(r==35)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[33]=0;f[33]<npar;f[33]++)
{      kcc=0;    for(k=0;k<33;k++) {  if(f[33]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[33]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[33]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[33] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[33]=0;g[33]<=1;g[33]++)

{   for(f[34]=0;f[34]<npar;f[34]++)
{      kcc=0;    for(k=0;k<34;k++)  {  if(f[34]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[34]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[34]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[34] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[34]=0;g[34]<=1;g[34]++)

{      if(r==37)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[35]=0;f[35]<npar;f[35]++)
{      kcc=0;    for(k=0;k<35;k++) {  if(f[35]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[35]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[35]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[35] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[35]=0;g[35]<=1;g[35]++)

{   for(f[36]=0;f[36]<npar;f[36]++)
{      kcc=0;    for(k=0;k<36;k++)  {  if(f[36]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[36]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[36]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[36] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[36]=0;g[36]<=1;g[36]++)

{      if(r==39)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[37]=0;f[37]<npar;f[37]++)
{      kcc=0;    for(k=0;k<37;k++) {  if(f[37]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[37]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[37]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[37] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[37]=0;g[37]<=1;g[37]++)

{   for(f[38]=0;f[38]<npar;f[38]++)
{      kcc=0;    for(k=0;k<38;k++)  {  if(f[38]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[38]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[38]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[38] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[38]=0;g[38]<=1;g[38]++)

{      if(r==41)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[39]=0;f[39]<npar;f[39]++)
{      kcc=0;    for(k=0;k<39;k++) {  if(f[39]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[39]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[39]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[39] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[39]=0;g[39]<=1;g[39]++)

{   for(f[40]=0;f[40]<npar;f[40]++)
{      kcc=0;    for(k=0;k<40;k++)  {  if(f[40]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[40]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[40]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[40] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[40]=0;g[40]<=1;g[40]++)

{      if(r==43)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[41]=0;f[41]<npar;f[41]++)
{      kcc=0;    for(k=0;k<41;k++) {  if(f[41]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[41]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[41]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[41] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[41]=0;g[41]<=1;g[41]++)

{   for(f[42]=0;f[42]<npar;f[42]++)
{      kcc=0;    for(k=0;k<42;k++)  {  if(f[42]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[42]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[42]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[42] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[42]=0;g[42]<=1;g[42]++)

{      if(r==45)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[43]=0;f[43]<npar;f[43]++)
{      kcc=0;    for(k=0;k<43;k++) {  if(f[43]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[43]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[43]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[43] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[43]=0;g[43]<=1;g[43]++)

{   for(f[44]=0;f[44]<npar;f[44]++)
{      kcc=0;    for(k=0;k<44;k++)  {  if(f[44]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[44]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[44]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[44] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[44]=0;g[44]<=1;g[44]++)

{      if(r==47)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[45]=0;f[45]<npar;f[45]++)
{      kcc=0;    for(k=0;k<45;k++) {  if(f[45]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[45]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[45]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[45] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[45]=0;g[45]<=1;g[45]++)

{   for(f[46]=0;f[46]<npar;f[46]++)
{      kcc=0;    for(k=0;k<46;k++)  {  if(f[46]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[46]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[46]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[46] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[46]=0;g[46]<=1;g[46]++)

{      if(r==49)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[47]=0;f[47]<npar;f[47]++)
{      kcc=0;    for(k=0;k<47;k++) {  if(f[47]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[47]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[47]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[47] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[47]=0;g[47]<=1;g[47]++)

{   for(f[48]=0;f[48]<npar;f[48]++)
{      kcc=0;    for(k=0;k<48;k++)  {  if(f[48]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[48]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[48]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[48] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[48]=0;g[48]<=1;g[48]++)

{      if(r==51)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[49]=0;f[49]<npar;f[49]++)
{      kcc=0;    for(k=0;k<49;k++) {  if(f[49]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[49]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[49]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[49] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[49]=0;g[49]<=1;g[49]++)

{   for(f[50]=0;f[50]<npar;f[50]++)
{      kcc=0;    for(k=0;k<50;k++)  {  if(f[50]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[50]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[50]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[50] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[50]=0;g[50]<=1;g[50]++)

{      if(r==53)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[51]=0;f[51]<npar;f[51]++)
{      kcc=0;    for(k=0;k<51;k++) {  if(f[51]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[51]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[51]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[51] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[51]=0;g[51]<=1;g[51]++)

{   for(f[52]=0;f[52]<npar;f[52]++)
{      kcc=0;    for(k=0;k<52;k++)  {  if(f[52]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[52]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[52]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[52] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[52]=0;g[52]<=1;g[52]++)

{      if(r==55)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[53]=0;f[53]<npar;f[53]++)
{      kcc=0;    for(k=0;k<53;k++) {  if(f[53]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[53]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[53]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[53] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[53]=0;g[53]<=1;g[53]++)

{   for(f[54]=0;f[54]<npar;f[54]++)
{      kcc=0;    for(k=0;k<54;k++)  {  if(f[54]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[54]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[54]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[54] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[54]=0;g[54]<=1;g[54]++)

{      if(r==57)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[55]=0;f[55]<npar;f[55]++)
{      kcc=0;    for(k=0;k<55;k++) {  if(f[55]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[55]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[55]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[55] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[55]=0;g[55]<=1;g[55]++)

{   for(f[56]=0;f[56]<npar;f[56]++)
{      kcc=0;    for(k=0;k<56;k++)  {  if(f[56]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[56]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[56]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[56] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[56]=0;g[56]<=1;g[56]++)

{      if(r==59)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[57]=0;f[57]<npar;f[57]++)
{      kcc=0;    for(k=0;k<57;k++) {  if(f[57]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[57]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[57]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[57] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[57]=0;g[57]<=1;g[57]++)

{   for(f[58]=0;f[58]<npar;f[58]++)
{      kcc=0;    for(k=0;k<58;k++)  {  if(f[58]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[58]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[58]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[58] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[58]=0;g[58]<=1;g[58]++)

{      if(r==61)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

    for(f[59]=0;f[59]<npar;f[59]++)
{      kcc=0;    for(k=0;k<59;k++) {  if(f[59]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[59]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[59]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++)   {     if(f[59] == cc[k])  {   kcc=1;   break;  }   }      if(kcc==1)   continue;
    for(g[59]=0;g[59]<=1;g[59]++)

{   for(f[60]=0;f[60]<npar;f[60]++)
{      kcc=0;    for(k=0;k<60;k++)  {  if(f[60]==f[k]) {   kcc=1;   break;  }   }    if(kcc==1)  continue;
       if(pares[f[60]][0]==qmag[(rg-r)/2][(rg+r-2)/2]||pares[f[60]][1]==qmag[(rg-r)/2][(rg+r-2)/2])  continue;
       for(k=0;k<=r-2; k++) {   if(f[60] == cc[k])  {   kcc=1;    break;  }   }     if(kcc==1)   continue;
    for(g[60]=0;g[60]<=1;g[60]++)

{      if(r==63)
       {  vp=0;
          colump=qmag[(rg+r-2)/2][(rg-r)/2]+pares[cc[0]][dd[0]]+2*ce-pares[f[r-3]][g[r-3]];
          for(k=0;k<r-3;k+=2)   colump+=2*ce-pares[f[k]][g[k]]+pares[f[k+1]][g[k+1]];
          if(colump!=r*qmag[(rg-1)/2][(rg-1)/2])  {   vp=1;    continue;   }
          op=0;     break;   }

   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }   if(op==0)  break;  }
   if(op==0)  break;  }   if(op==0)  break;  }
   }

