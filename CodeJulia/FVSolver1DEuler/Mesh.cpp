#include "Mesh.h"

Mesh::Mesh( int nCells, Settings* settings ) : _nCells( nCells ), _settings( settings ) {}

Mesh::~Mesh() {}

int Mesh::GetNumberOfCells() const { return _nCells; }
