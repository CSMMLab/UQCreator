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
        double a = settings->get_as<double>("a").value_or(0.0);
        double b = settings->get_as<double>("b").value_or(0.0);
        if((a==0 && b==0) || _dimension == -1 || _numCells == -1)
            std::cerr << "[ERROR]: Mesh::Mesh invalid mesh parameters" << std::endl;
        else{
            this->CreateGrid(a,b);
            _status = MESH_STATUS_LOADED;
        }
    }
}

Mesh::~Mesh(){

}

void Mesh::Load(std::string filename){
    std::cerr << "[ERROR]: Mesh::Load not yet implemented" << std::endl;
}

void Mesh::CreateGrid(double a, double b){
    assert(_dimension == 1);
    _mesh.resize(_numCells+4);
    _spacing.resize(_numCells+4-1);
    for(int j=-2; j<_numCells+2; ++j){
        _mesh[j+2] = a + j*(b-a)/(_numCells-1);
        if(j!=0){
            _spacing[j+1] = _mesh[j+2] - _mesh[j+1];
        }
    }
}

int Mesh::GetNumCells() const{
    return _numCells;
}

int Mesh::GetDimension() const{
    return _dimension;
}

blaze::DynamicVector<double> Mesh::GetGrid(){
    return _mesh;
}

blaze::DynamicVector<double> Mesh::GetSpacing(){
    return _spacing;
}

