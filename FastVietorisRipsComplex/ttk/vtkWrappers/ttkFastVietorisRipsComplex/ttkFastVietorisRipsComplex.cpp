#include                  <ttkFastVietorisRipsComplex.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <chrono>
#include <thread>
#include <math.h>
#include <cmath>
#include <vector>

vtkStandardNewMacro(ttkFastVietorisRipsComplex)
//struct edge
//{
//    int data;
//    edge *nextEdge;
//};
//struct node
//{
//    int data;
//    node *next;
//    edge *edge;
//};
////class list
////{
////private:
////    node *head, *tail;
////public:
////    list()
////    {
////        head=NULL;
////        tail=NULL;
////    }
////};
// node *head, *tail;
int ttkFastVietorisRipsComplex::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
    vtkPointSet *input =vtkPointSet::SafeDownCast (inputs[0]);
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
    vtkPoints *inputPointSet = (input)->GetPoints();
    edges.clear();
    tetras.clear();
    triangles.clear();

    
    int size = input->GetNumberOfPoints();
    vector<vector<int> > neighbors;
    neighbors.reserve(size);
    Timer t;
    int k=0;
    if (Option3) {
        k=3;
    }
    else if (Option2) {
        k=2;
    }
    else if (Option1) {
        k=1;
    }
    cout<<"k= "<<k<<endl;

    for (int i=0; i<size; i++) {
        vector<int> graph;
        for (int j=0; j<i; j++) {
            double * mypointer;
            mypointer= input->GetPoint(i);
            double a[3];
            a[0]=mypointer[0];
            a[1]=mypointer[1];
            a[2]=mypointer[2];
            double * mypointer2;
            mypointer2= input->GetPoint(j);
            double a2[3];
            a2[0]=mypointer2[0];
            a2[1]=mypointer2[1];
            a2[2]=mypointer2[2];
            double distance = 0;

            for(int i = 0; i < 3; i++){
                distance += (a[i] - a2[i])*(a[i] - a2[i]);
            }
            distance=sqrt(distance);
            if (distance <= 2*Radius) {

                graph.push_back(j);
            }
        }
        neighbors.push_back(graph);
    }
//    for (int i=0; i<size; i++) {
//        node *graph=new node;
//        graph->data=-1;
//        graph->next=head;
//        head=graph;
//        for (int j=0; j<i; j++) {
//            double * mypointer;
//            mypointer= input->GetPoint(i);
//            double a[3];
//            a[0]=mypointer[0];
//            a[1]=mypointer[1];
//            a[2]=mypointer[2];
//            double * mypointer2;
//            mypointer2= input->GetPoint(j);
//            double a2[3];
//            a2[0]=mypointer2[0];
//            a2[1]=mypointer2[1];
//            a2[2]=mypointer2[2];
//            double distance = 0;
//
//            for(int i = 0; i < 3; i++){
//                distance += (a[i] - a2[i])*(a[i] - a2[i]);
//            }
//            distance=sqrt(distance);
//            if (distance <= 2*Radius) {
//                edge *e=new edge;
//                e->data=j;
//                e->nextEdge=graph->edge;
//                graph->edge=e;
//                //graph.push_back(j);
//            }
//        }
//        //neighbors.push_back(graph);
//    }
    
    {
        stringstream msg;
        msg << "[FastVietorisRipsComplex] Graph processed in "<< t.getElapsedTime() << " s."<< endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Timer t2;
    incrementalVR(neighbors,k);
    {
        stringstream msg;
        msg << "[FastVietorisRipsComplex] VR processed in "<< t2.getElapsedTime() << " s."<< endl;
        dMsg(cout, msg.str(), timeMsg);
    }
    {
        stringstream msg;
        msg << "[FastVietorisRipsComplex] FastVietorisRipsComplex done in "<< t.getElapsedTime() << " s."<< endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    vtkUnstructuredGrid * grid=vtkUnstructuredGrid::New();
    grid->Allocate();
    grid->SetPoints(inputPointSet);
    
    if(Option1){
        cout<<"EDGES "<<edges.size()<<endl;
        vtkIdType ids[2];
        for (int i=0; i<edges.size(); i++) {
            ids[0]=(edges.at(i)[0]);
            ids[1]=(edges.at(i)[1]);
            grid->InsertNextCell(VTK_LINE,2,ids);

        }
    }
    if(Option2){
        cout<<"TRIANGLES "<<triangles.size()<<endl;
        vtkIdType Triangleids[3];
        for (int i=0; i<triangles.size(); i++) {
            Triangleids[0]=(triangles.at(i)[0]);
            Triangleids[1]=(triangles.at(i)[1]);
            Triangleids[2]=(triangles.at(i)[2]);
            grid->InsertNextCell(VTK_TRIANGLE,3,Triangleids);

        }
    }
    
    if(Option3){
        cout<<"TETRAS "<<tetras.size()<<endl;
        vtkIdType TetraIds[4];
        for (int i=0; i<tetras.size(); i++) {
            TetraIds[0]=(tetras.at(i)[0]);
            TetraIds[1]=(tetras.at(i)[1]);
            TetraIds[2]=(tetras.at(i)[2]);
            TetraIds[3]=(tetras.at(i)[3]);
            grid->InsertNextCell(VTK_TETRA,4,TetraIds);

        }
    }
    output->ShallowCopy(grid);
    
    
    
  {
    stringstream msg;
    msg << "[ttkFastVietorisRipsComplex] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}

void ttkFastVietorisRipsComplex::incrementalVR(const vector<vector<int> > &neighbors,int k){
    for (int i = 0; i < neighbors.size(); i++) {
        vector<int> index;
        index.push_back(i);
        addCofaces(neighbors,k, index, neighbors[i]);
    }
}

void ttkFastVietorisRipsComplex::addCofaces(const vector<vector<int> > &neighbors,int k, const vector<int> &index,const vector<int> &local){
    
    if (index.size()==2) {
        edges.push_back(index);
    }
    else if (index.size()==3) {
        triangles.push_back(index);
    }
    else if (index.size()==4) {
        tetras.push_back(index);
    }
    if (index.size() > k) {
        return;
    }
    else{
        for (int i = 0; i < local.size(); i++) {
            int v = local[i];
            vector<int> s;
            s.push_back(v);
            for (int l=index.size()-1; l>=0; l--) {
                s.push_back(index[l]);
            }
            vector<int> M;
            intersection(local, neighbors[v],M);
            addCofaces(neighbors,k, s, M);
        }
    }
}

void ttkFastVietorisRipsComplex::intersection(const vector<int> &a,const vector<int> &b,vector<int> &result){
    
    int ai = 0;
    int bi = 0;
    
    
    while (ai < a.size() && bi < b.size()) {
        if (a[ai] < b[bi]) {
            ai++;
        } else if (a[ai] > b[bi]) {
            bi++;
        } else {
            result.push_back(a[ai]);
            ai++;
            bi++;
        }
    }
    
}

 

int ttkFastVietorisRipsComplex::FillOutputPortInformation(int port, vtkInformation* info){
    if(port == 0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    
    return 1;
}
int ttkFastVietorisRipsComplex::FillInputPortInformation(int port, vtkInformation* info){
    if(port == 0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
    
    return 1;
}
