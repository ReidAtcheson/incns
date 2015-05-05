#Very simple visualization output.
#Uses template from http://www.paraview.org/Wiki/ParaView/Data_formats#Displaying_data_as_structured_grid
function tocsv(xs,ys,zs,fn,filename)
    f=open(filename,"w");
    print(f,"x coord,y coord,z coord,field\n");
    for x in xs
        for y in ys
            for z in zs
                of=map(abs,fn(x,y,z));
                out=hypot(hypot(of[1],of[2]),of[3]);
                print(f,"$(x), $(y), $(z), $(out)\n");
            end
        end
    end
    close(f);
end
