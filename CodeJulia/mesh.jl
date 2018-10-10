__precompile__

struct Mesh
    # number of cells
    Nx::Int64;
    Ny::Int64;

    # Array of mesh nodes positions
    nodes::Array{Float64,3};
    # Array of cell nodes
    cellNodes::Array{Float64,4};
    # Array of cell midpoints
    cellMidPoints::Array{Float64,3};
    # Array of cell Volumes
    cellVolumes::Array{Float64,2};
    # Array of face normals
    faceNormals::Array{Float64,4};
    # Array of face normals
    faceUnitNormals::Array{Float64,4};

    function Mesh(settings::Settings)
        cellType::Int = 4; # number of faces per cell
        dimension::Int = 2;
        x = linspace(settings.a,settings.b,settings.Nx+1);
        y = linspace(settings.c,settings.d,settings.Ny+1);
        dx = x[2]-x[1];
        dy = y[2]-y[1];
        nodes = x*y';
        nodes = zeros(settings.Nx+1,settings.Ny+1,dimension);
        cellNodes = zeros(settings.Nx,settings.Ny,cellType,dimension);
        cellMidPoints = zeros(settings.Nx,settings.Ny,dimension);
        cellVolumes = zeros(settings.Nx,settings.Ny);
        faceNormals = zeros(settings.Nx,settings.Ny,cellType,dimension);
        faceUnitNormals = zeros(settings.Nx,settings.Ny,cellType,dimension);

        # write nodes
        for i = 1:(settings.Nx+1)
            for j = 1:(settings.Ny+1)
                nodes[i,j,1] = x[i];
                nodes[i,j,2] = y[j];
            end
        end
        # write cell properties
        for i = 1:settings.Nx
            for j = 1:settings.Ny                           #        s_2
                cellNodes[i,j,1,:] = nodes[i,j,:];          #       2------3
                cellNodes[i,j,2,:] = nodes[i,j+1,:];        #       |      |  s_3
                cellNodes[i,j,3,:] = nodes[i+1,j+1,:];      #   s_1 |      |  cell (i,j)
                cellNodes[i,j,4,:] = nodes[i+1,j,:];        #       1------4
                                                            #         s_4
                cellVolumes[i,j] = 0.5*( (nodes[i+1,j+1,1]-nodes[i,j,1])*(nodes[i,j+1,2]-nodes[i+1,j,2])-
                                        (nodes[i+1,j+1,2]-nodes[i,j,2])*(nodes[i,j+1,1]-nodes[i+1,j,1]) );

                cellMidPoints[i,j,:] = 0.25*(cellNodes[i,j,1,:]+cellNodes[i,j,2,:]+cellNodes[i,j,3,:]+cellNodes[i,j,4,:]);

                # compute normals for each face
                for l = 1:cellType
                    t = cellNodes[i,j,l%cellType+1,:]-cellNodes[i,j,l,:];
                    # make sure all faceNormals point into the positive direction
                    if t[2] >= 0.0 && t[1] <= 0.0 
                        n = [t[2];-t[1]];
                    else
                        n = [-t[2];t[1]];
                    end
                    faceNormals[i,j,l,:] = n;
                    faceUnitNormals[i,j,l,:] = n/norm(n);
                    #println("dy = ",dy,", dx = ",dx," ,norm of faceNormal is ",norm(faceNormals[i,j,:]));
                end

            end
        end

        new(settings.Nx,settings.Ny,nodes,cellNodes,cellMidPoints,cellVolumes,faceNormals,faceUnitNormals);
    end

end

function InsertObstacle(obj::Mesh, xM, radius, flowVolume)
    #flowVolume = zeros(obj.Nx,obj.Ny)
    for i = 1:obj.Nx
        for j = 1:obj.Ny
            if norm(obj.cellMidPoints[i,j,:]-xM) < radius
                flowVolume[i,j] = 1;
            end
        end
    end
    return flowVolume;
end

function InsertRectangularObstacle(obj::Mesh, xM, radius, flowVolume)
    #flowVolume = zeros(obj.Nx,obj.Ny)
    for i = 1:obj.Nx
        for j = 1:obj.Ny
            if max(abs(obj.cellMidPoints[i,j,1]-xM[1]),abs(obj.cellMidPoints[i,j,2]-xM[2])) < radius
                flowVolume[i,j] = 1;
            end
        end
    end
    return flowVolume;
end

function InsertRectangularObstacle(obj::Mesh, xM, lengthX,lengthY, flowVolume)
    for i = 1:obj.Nx
        for j = 1:obj.Ny
            if abs(obj.cellMidPoints[i,j,1]-xM[1]) < lengthX && abs(obj.cellMidPoints[i,j,2]-xM[2]) < lengthY
                flowVolume[i,j] = 1;
            end
        end
    end
    return flowVolume;
end

function AddWallBoundary(obj::Mesh,wall::String, flowVolume)
    #flowVolume = zeros(obj.Nx,obj.Ny)

    if wall == "north" || wall == "south"
        for i = 1:obj.Nx
            if wall == "north"
                flowVolume[i,end] = 1;
            else
                flowVolume[i,1] = 1;
            end
        end
    end
    if wall == "west" || wall == "east"
        for i = 1:obj.Ny
            if wall == "west"
                flowVolume[1,i] = 1;
            else
                flowVolume[end,i] = 1;
            end
        end
    end
    return flowVolume;
end