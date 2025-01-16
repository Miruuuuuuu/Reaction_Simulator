#include<iostream>
using namespace std;
void merge(int arr[],int low,int high)
{
   int s=(high-low)+1;
   int m=(high+low)/2;
   int k=0;
   int i=low;
   int j=m+1;
   int *temp=new int[s];
while(i<=m && j<high)
{
    if(arr[i]<arr[j])

    {
       temp[k++]=arr[i++]; 
    }
    else if(arr[i]>arr[j])
    {
        temp[k++]=arr[i++];
    }
}
while(i<=m)
{
    temp[k++]=arr[i++];
}
while (j<high)
{
    temp[k++]=arr[j++];
}
for(int i=0;i<s;i++)
{
    arr[low+1]=temp[i];
}
delete[] temp;
}
void subarray(int arr[],int s,int e)
{
    if (s<e)
    {
        int m=(s+e)/2;
        subarray(arr,s,m);
        subarray(arr,m+1,e);
        merge(arr, s, e);   
     }
    
    
    
}
int main()
{
int arr[5]={21,22,1,2,3};
int s=0;
int e=4;

subarray(arr,s,e);
for(int i=0;i<5;i++)
{
    cout<<arr[i]<<endl;
}

}