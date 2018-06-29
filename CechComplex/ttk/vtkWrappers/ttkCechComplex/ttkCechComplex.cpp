#include                  <ttkCechComplex.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <chrono>
#include <thread>
#include <math.h>
#include <cmath>
#include <vector>

vtkStandardNewMacro(ttkCechComplex)

int ttkCechComplex::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
    
   input =vtkPointSet::SafeDownCast (inputs[0]);
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
//
    Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;

  triangulation->setWrapper(this);
  cechComplex_.setupTriangulation(triangulation);
  cechComplex_.setWrapper(this);
    edges.clear();
    tetras.clear();
    triangles.clear();
    
    int size1 = input->GetNumberOfPoints();
    vector<int> v2(1, -1);
    vector<vector<int> > neighbors(size1,v2);
    
   // neighbors(input->GetNumberOfPoints());
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  vtkPoints *inputPointSet = (input)->GetPoints();
    
    for (int i=0; i<size1-1; i++) {
        for (int j=i+1; j<size1; j++) {
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
            if (intersect(distance,Radius)) {

                vector<int> edge;
                edge.push_back(i);
                edge.push_back(j);
                edges.push_back(edge);
                //cout<<i<<" "<<j<<endl;
                if (neighbors[i][0] == -1) {
                    neighbors[i][0]=j;
                }
                else{
                    neighbors[i].push_back(j);}
            }
            
        }
        
    }
    
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
    KSimplices(neighbors,k);
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
    msg << "[ttkCechComplex] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
int ttkCechComplex::FillOutputPortInformation(int port, vtkInformation* info){
    if(port == 0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    
    return 1;
}
int ttkCechComplex::FillInputPortInformation(int port, vtkInformation* info){
    if(port == 0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
    
    return 1;
}
bool ttkCechComplex::intersect(double distance,double Radius){
    // Circles do not overlap
    if (distance > 2*Radius) {
        return false;
    }
    // One circle contains the other
    if (distance < 0) {
        return false;
    }
     //These are the same circle
    if (distance == 0) {
        return false;
    }

    return true;
}
void ttkCechComplex::KSimplices(const vector<vector<int> > &neighbors,int max){
    int k = 2;
    vector<int> candidates;
    vector<vector<int>> Sk;

    while (1) {
        Sk.clear();
        for (int i = 0; i < input->GetNumberOfPoints(); i++) {
            combination.clear();
            if (neighbors[i][0]!= -1) {
                    getCombination(neighbors[i], i, k);
            }
            
            for (int j = 0; j < combination.size(); j++) {
                vector<int> u = combination[j];
                if (verify(u)) {
                    Sk.push_back(u);
                }
            }
//            cout<<"==== "<<endl;
//            std::chrono::seconds dura(10);
//            std::this_thread::sleep_for( dura );
//
            
        }
        if (Sk.size()!=0) {
            for (int z=0; z<Sk.size(); z++) {
                //S.push_back(Sk[z]);
                if (Sk[z].size()==3) {
                    triangles.push_back(Sk[z]);
                }
                else if (Sk[z].size()==4){
                        tetras.push_back(Sk[z]);}
            }
            k++;
        } else {
            break;
        }
//
        if (k > max) {
            break;
        }



    }


}



bool ttkCechComplex::verify(const vector<int> &arr){
    
    int k = arr.size();
    vector<Point>uCells;

    for (int i = 0; i < k-1; i++) {
        for (int j = 1; j < k; j++) {
            vector<Point> item= getIntersection(i,j,arr);
            for (int l=0; l<item.size(); l++) {
                if (item[i].x!=-1) {
                    uCells.push_back(item[i]);
                }
            }
            
        }
    }
    int i= uCells.size();

    while (i--) {
        Point p = uCells[i];
        bool exists = true;
        int j = arr.size();
        while (j--) {
            double * mypointer;
            mypointer= input->GetPoint(arr[j]);
            Point cell;
            cell.x=mypointer[0];
            cell.y=mypointer[1];
            cell.z=mypointer[2];
            if (!isPointInsideCell(p, cell)) {
                exists = false;
                break;
            }

        }
        if (exists) {
            return true;
        }
    }
    return false;
    
}


bool ttkCechComplex::isPointInsideCell(Point point,Point cell){
    double distance = 0;
    distance += (point.x - cell.x)*(point.x - cell.x);
    distance += (point.y - cell.y)*(point.y - cell.y);

    distance=sqrt(distance);
    if(distance<= Radius)
        return true;
    else
        return false;
    
    
    
}


vector<Point> ttkCechComplex::getIntersection(int i,int j,const vector<int> &arr){
    
        double * mypointer;
        vector<Point> result;
        Point p;
        Point p2;

        mypointer= input->GetPoint(arr[i]);
        double a[3];
        a[0]=mypointer[0];
        a[1]=mypointer[1];
        a[2]=mypointer[2];
    
        double * mypointer2;
        mypointer2= input->GetPoint(arr[j]);
        double a2[3];
        a2[0]=mypointer2[0];
        a2[1]=mypointer2[1];
        a2[2]=mypointer2[2];
    
        double distance = 0;
        for(int i = 0; i < 3; i++){
            distance += (a[i] - a2[i])*(a[i] - a2[i]);
        }
    
        distance=sqrt(distance);
        if (intersect(distance,Radius)==false) {
            p.x=-1;
            result.push_back(p);
            result.push_back(p);
            return result;
        }
        double aa = (pow(Radius, 2) - pow(Radius, 2) + pow(distance, 2)) / (2 * distance);
        double h = sqrt(pow(Radius, 2) - pow(aa, 2));
        double temp[2] = {
                     a[0] + aa * (a2[0] - a[0]) / distance,
                     a[1] + aa * (a2[1] - a[1]) / distance//,
                    // a[2] + aa * (a2[2] - a[2]) / distance
        };
    
        p.x=temp[0] + h * (a2[1] - a[1]) / distance;
        p.y=temp[1] - h * (a2[0] - a[0]) / distance;
    
        p2.x=temp[0] - h * (a2[1] - a[1]) / distance;
        p2.y=temp[1] + h * (a2[0] - a[0]) / distance;
    
        result.push_back(p);
        result.push_back(p2);

    
        return result;

}
void ttkCechComplex::getCombination(const vector<int> &arr,int index, int r){

    int n = arr.size();
    setCombination(arr, n, r);
    for (int i = 0; i < combination.size(); i++) {
        combination[i].insert(combination[i].begin(),index);
    }
}
void ttkCechComplex::setCombination(const vector<int> &arr, int n, int r){
    vector<int> sol(r);
    allCombinations(arr, sol, 0, n-1, 0, r);
}

void ttkCechComplex::allCombinations(const vector<int> &arr, vector<int> &sol, int start, int end,int index, int r){
    if (index == r)
    {
        combination.push_back(sol);
        return;
    }

    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        sol[index] = arr[i];
        allCombinations(arr, sol, i+1, end, index+1, r);
    }
}
