#p_era=fERA['lsp'][:,:,:]
#x_era=fERA['longitude'][:]
#y_era=fERA['latitude'][:]

function gridPrecip(pmap_era,cmap_era,pmap_imerg,cmap_imerg,x1,y1,p_era,
                    p_imerg,it)
    np=size(x1)[2]
    #println(size(x1))
    p_era1d=[]
    for itime=1:40
        for ip=1:np
            x11=x1[itime,ip]
            if x11<0
                x11=x11+360
            end
            y11=y1[itime,ip]
            iy=Int(trunc((90-y11)/0.75))
            ix=Int(trunc(x11/0.75))
            fy=(90-y11)/0.75-iy
            fx=x11/0.75-ix
            ix=ix+1
            iy=iy+1
            ix1=ix+1
            iy1=iy+1
            if ix1==481
                ix1=1
            end
            if ix==481
                ix=1
                ix1=2
            end
            p=-999
            try
                p=(1-fx)*(1-fy)*p_era[itime+it,iy,ix]+
                (1-fx)*(fy)*p_era[itime+it,iy,ix1]+
                (fx)*(1-fy)*p_era[itime+it,iy1,ix]+
                (fx)*(fy)*p_era[itime+it,iy1,ix1]
            catch
                println("$(ix) $(iy)")
                exit(1)
            end
            try
                ix=Int(trunc(x11))+1
                if ix==361
                    ix=1
                end
                iy=Int(trunc(y11+90))+1
                if iy>180
                    iy=180
                end
                #if p>0
                pmap_era[ix,iy]=pmap_era[ix,iy]+1000*p
                cmap_era[ix,iy]=cmap_era[ix,iy]+1
                #end
            catch
                println("gridding $(x11) $(y11) $(ix) $(iy)")
                exit(1)
            end
            push!(p_era1d,p)
        end
    end
    p_imerg1d=[]
    for itime=1:40
        for ip=1:np
            x11=x1[itime,ip]
            y11=y1[itime,ip]
            iy=Int(trunc((y11+89.7005)/0.6))
            ix=Int(trunc((x11+179.7)/0.6))
            fy=(y11+89.7005)/0.6-iy
            fx=((x11+179.7)/0.6)-ix
            ix=ix+1
            iy=iy+1
            ix1=ix+1
            iy1=iy+1
            if ix1==601
                ix1=1
            end
            if iy==300
                iy1=300
            end
            p=-999.9
            try
                p=(1-fx)*(1-fy)*p_imerg[itime+it,ix,iy]+
                (1-fx)*(fy)*p_imerg[itime+it,ix,iy1]+
                (fx)*(1-fy)*p_imerg[itime+it,ix1,iy]+
                (fx)*(fy)*p_imerg[itime+it,ix1,iy1]
            catch
                println("$(fx) $(fy) $(ix) $(iy) $(y11)")
            end
            if x11<0
                x11=x11+360
            end
            ix=Int(trunc(x11))+1
            iy=Int(trunc(y11+90))+1
            if ix==361
                ix=1
            end
            iy=Int(trunc(y11+90))+1
            if iy>180
                iy=180
            end
            if p>-0.1 
                pmap_imerg[ix,iy]=pmap_imerg[ix,iy]+p
                cmap_imerg[ix,iy]=cmap_imerg[ix,iy]+1
            end
            push!(p_imerg1d,p)

        end
      end
    
    println(sum(cmap_era))
    return p_era1d,p_imerg1d,cmap_era,pmap_era,
    cmap_imerg,pmap_imerg
end
