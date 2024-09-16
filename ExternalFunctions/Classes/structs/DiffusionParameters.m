classdef DiffusionParameters
    properties
        dx
        nPts
        filter
        filterSize
        flag3d
        dz
        expTime
        locPrecXY
        locPrecZ
        frametime
    end
    methods
        function obj=DiffusionParameters(dx, numpts, filter, filtersi, flag3d, dz, expT, locPXY, locPZ, frametime)
            obj.dx = dx;
            obj.nPts = numpts;
            obj.filter = filter;
            obj.filterSize = filtersi;
            obj.flag3d = flag3d;
            obj.dz = dz;
            obj.expTime = expT;
            obj.locPrecXY = locPXY;
            obj.locPrecZ = locPZ;
            obj.frametime = frametime;
        end
        
        function res = isequal(object)
            if obj == object
                res = 1;
            else
                res = 0;
            end
        end
    end
end