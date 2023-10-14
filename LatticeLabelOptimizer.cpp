// C++ program to distinct permutations of the string 
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "tensor_type.h"
#include "Matrix.h"
using namespace std;


int g_value( int Length_, Mat_2_int C_){

int temp=0;

for(int i=0;i<Length_;i++){
for(int j=0;j<Length_;j++){
temp += C_[i][j]*((j-i))*(j-i);//*abs(j-i)*abs(j-i);
}}

return temp;
}

int N_for_given_partition(int l, int Length_, Mat_2_int C_){
int temp=0;
for(int i=0;i<=l;i++){
for(int j=l+1;j<=Length_-1;j++){
temp += C_[i][j];
}}
return temp;
}



int D_for_given_partition(int l, int Length_, Mat_2_int C_){
int Dl_max=-1000;
int temp;

for(int i=0;i<=l;i++){
for(int j=l+1;j<=Length_-1;j++){
temp = (C_[i][j])*abs(j-i);
if(temp>Dl_max){
Dl_max=temp;
}
}}
return Dl_max;
}


int Xc_for_the_givenC(int Length_, Mat_2_int C_, string TYPE){
int Xl_max=-1000;
int temp;
for(int l=0;l<=Length_-2;l++){
if(TYPE=="N"){
temp=N_for_given_partition(l, Length_, C_);
}
else if(TYPE=="D"){
temp=D_for_given_partition(l, Length_, C_);
}
else{
cout<<"Only TYPE D or N is allowed"<<endl;
assert(false);
}
if(temp>Xl_max){
Xl_max=temp;}
}
return Xl_max;
}




int main(){

int lx=8;
int ly=10;
int Length=100;//lx*ly;
Mat_2_int C;

double Temperature;
double beta;

C.resize(Length);
for(int i=0;i<Length;i++){
C[i].resize(Length);
for(int j=0;j<Length;j++){
C[i][j]=0;
}
}


//For PBC chain

for(int i=0;i<=Length-2;i++){
C[i][i+1]=1;
C[i+1][i]=1;
}
C[0][Length-1]=1;
C[Length-1][0]=1;




//For 2d OBC
/*
int i_,neigh;
for(int ix=0;ix<lx;ix++){
for(int iy=0;iy<ly;iy++){
i_= ix+lx*iy;
//+x
neigh = (ix+1)+lx*iy;
if(ix < lx-1){
C[i_][neigh]=1;
C[neigh][i_]=1;
}
//+y
neigh = (ix) + lx*(iy+1);
if(iy < ly-1){
C[i_][neigh]=1;
C[neigh][i_]=1;
}
}}
*/

//For 2X5 ladder, rung coupled
/*
for(int i=0;i<(Length/2);i++){
C[i][i+(Length/2)]=1;
C[i+(Length/2)][i]=1;
}
*/
/*
for(int i=0;i<(Length/2);i++){
C[i][Length-1-i]=1;
C[Length-1-i][i]=1;
}
*/

cout<<"Initial Nc = "<<Xc_for_the_givenC(Length, C, "N")<<endl;
cout<<"Initial Dc = "<<Xc_for_the_givenC(Length, C, "D")<<endl;
cout<<"Initial g_value = "<<g_value(Length, C)<<endl;


cout<<"Printing C:"<<endl;
for(int i=0;i<Length;i++){
for(int j=0;j<Length;j++){
cout<<C[i][j]<<" ";
}
cout<<endl;
}
cout<<"--------------"<<endl;


Mat_1_int P_inv,P;
P_inv.resize(Length);
P.resize(Length);
for(int i=0;i<Length;i++){
P_inv[i]=i;
P[i]=i;
}

Mat_2_int Cp_inv, Cp_inv_old;
Cp_inv.resize(Length);
Cp_inv_old.resize(Length);
for(int i=0;i<Length;i++){
Cp_inv[i].resize(Length);
Cp_inv_old[i].resize(Length);
for(int j=0;j<Length;j++){
Cp_inv[i][j]= C[ P_inv[i]  ][ P_inv[j] ];
Cp_inv_old[i][j]= C[ P_inv[i]  ][ P_inv[j] ];
}
}

int seed=46;
srand(seed);

int Nc_old, Nc_, Dc_old, Dc_, g_old, g_;


Nc_old = Xc_for_the_givenC(Length, Cp_inv, "N");
Dc_old = Xc_for_the_givenC(Length, Cp_inv, "D");
g_old = g_value(Length, Cp_inv);

//Mat_1_doub Temperature_array = {1000.0, 800.0, 600.0, 400.0, 200.0, 100.0, 70.0, 50.0, 30.0, 20.0, 15.0, 14.0, 13.0, 12.0, 10.0, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.1, 0.01};

Mat_1_doub Temperature_array = {500.0, 400.0, 300.0, 200.0, 100.0, 70.0, 50.0, 30.0, 20.0, 15.0, 14.0, 12.0, 10.0, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.75, 5.5, 5.25, 5.0, 4.75, 4.5, 4.25, 4.0, 3.75, 3.5, 3.25, 3.0, 2.75, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.5, 0.1, 0.01};


for(int Temperature_ind=0;Temperature_ind<52;Temperature_ind++){
Temperature = Temperature_array[Temperature_ind];
beta = 1.0/Temperature;

//--------------- Optimizer---------------------//


int Max_Flips_tried=100000;
int flips_done=0;
int flips_tried=0;

int temp2, temp1;
bool flip_accepted;
bool flip_is_junk;
int label1, label2, label1_old, label2_old;
int c1, c2, c1_old, c2_old; //coordinates
while(flips_tried<=Max_Flips_tried){

/* generate 2 integers b/w 0 and L-1 */
int c1 = rand() % (Length);
int c2 = rand() % (Length);

if( ((c1 == c1_old) && (c2 == c2_old))  ||
    ((c1 == c2_old) && (c2 == c1_old))	||
    (c1==c2)
   ){
flip_is_junk=true;
}
else{
flip_is_junk=false;
}



if(!flip_is_junk){

cout<<c1<<"  "<<c2<<endl;

label1=P[c1];
label2=P[c2];

assert(label1!=label2);

for(int i=0;i<Length;i++){
Cp_inv[i][label1] = Cp_inv_old[i][label2];
Cp_inv[i][label2] = Cp_inv_old[i][label1];
}
for(int i=0;i<Length;i++){
temp2=Cp_inv[label2][i];
temp1=Cp_inv[label1][i];
Cp_inv[label1][i] = temp2;
Cp_inv[label2][i] = temp1;
}

Nc_ = Xc_for_the_givenC(Length, Cp_inv, "N");
Dc_ = Xc_for_the_givenC(Length, Cp_inv, "D");
g_ = g_value(Length, Cp_inv);

double p = abs(rand()/(1.0*RAND_MAX));
//if(Nc_<Nc_old){
//if(g_<g_old){
//if( (Dc_<Dc_old && Nc_==Nc_old)){

//cout<<"check : "<<exp(1.0*beta*(1.0*(g_-g_old)))<<"   "<<p<<endl;
//cout<<"check2 : "<<beta<<"  "<<g_<<"  "<<g_old<<endl;

if(exp(-1.0*beta*(g_-g_old))>p){// && Nc_<=Nc_old){

//flip is accepted;
P[c1]=label2;
P[c2]=label1;
P_inv[label1]=c2;
P_inv[label2]=c1;
flips_done++;
Nc_old=Nc_;
Dc_old=Dc_;
Cp_inv_old=Cp_inv;
g_old=g_;
}
else{
//flip is not accepted;
//Dont change Nc_old and Cp_inv_old
Nc_=Nc_old; //just to write the Nc_ for present configuration
Dc_=Dc_old;
P[c1]=label1;
P[c2]=label2;
P_inv[label1]=c1;
P_inv[label2]=c2;
Cp_inv = Cp_inv_old;
g_=g_old;
}

cout<<"Stats : "<<flips_tried<<"  "<<flips_done<<"  "<<Nc_<<"  "<<Dc_<<"  "<<g_<<endl;
cout<<"P : ";
for(int i=0;i<Length;i++){
cout<<P[i]<<" ";
}
cout<<endl;
/*
cout<<"Printing C:"<<endl;
for(int i=0;i<Length;i++){
for(int j=0;j<Length;j++){
cout<<Cp_inv[i][j]<<" ";
}
cout<<endl;
}
*/
cout<<"----------"<<endl;




flips_tried++;
}

}



cout<<"Final Nc(T="<< Temperature <<") = "<<Xc_for_the_givenC(Length, Cp_inv, "N")<<endl;
cout<<"Final Dc(T="<< Temperature <<") = "<<Xc_for_the_givenC(Length, Cp_inv, "D")<<endl;
cout<<"Final g(T="<< Temperature <<") = "<<g_value(Length, Cp_inv)<<endl;

cout<<"Printing C for Temperature = : "<< Temperature<<endl;
for(int i=0;i<Length;i++){
for(int j=0;j<Length;j++){
cout<<Cp_inv[i][j]<<" ";
}
cout<<endl;
}


}//Temperature

return 0;
}
