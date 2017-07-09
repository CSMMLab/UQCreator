#include "mesh.h"

Mesh::Mesh(std::string inputFile) : _status(MESH_STATUS_UNLOADED){
    auto file = cpptoml::parse_file(inputFile);
    auto settings = file->get_table("mesh");
    std::string meshFile = settings->get_as<std::string>("fileName").value_or("none");
    if(meshFile.compare("none")){
        this->Load(meshFile);
    }
    else{
        _dimension = settings->get_as<int>("Dimension").value_or(-1);
        _numCells = settings->get_as<int>("NumberOfCells").value_or(-1);
        _numCells = settings->get_as<int>("NumberOfCells").value_or(-1);
        double a = settings->get_as<int>("b").value_or(0);
        double b = settings->get_as<int>("b").value_or(0);
        if(a!=0 && b!=0)
            this->CreateGrid(a,b);
    }
}

void Mesh::Load(std::string filename){

}

void Mesh::CreateGrid(double a, double b){
    _mesh.resize(_dimension);
    _spacing.resize(_dimension);
    for(int i=0; i<_dimension; ++i){
        _mesh[i] = vector(_numCells);
        _spacing[i] = vector(_numCells-1);
    }
}

int Mesh::GetNumCells() const{
    return _numCells;
}

int Mesh::GetDimension() const{
    return _dimension;
}

meshData& Mesh::GetGrid(){
    return _mesh;
}

meshData& Mesh::GetSpacing(){
    return _spacing;
}

