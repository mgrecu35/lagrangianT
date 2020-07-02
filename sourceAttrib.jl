function attrib(p,rho,x,y,T,q,dp,src,csrc,old2d,cold2d)
    #src=zeros(360,90,20)
    #csrc=zeros(360,90,20)
    #oldM=zeros(360,90,20)
    println(size(p))
    println(size(src))
    nt,np=size(p)
    println(nt)
    oldM=0.0
    evapM=0.0
    for i=1:np
        fract=Array{Float32}(undef,0)
        ix=Array{Float32}(undef,0)
        iy=Array{Float32}(undef,0)
        ip=Array{Float32}(undef,0)
        itime=Array{Int32}(undef,0)
        qs=Array{Float32}(undef,0)
        old=1.0
        qold=q[nt,i]
        for j=nt:-1:2
            if q[j,i]*rho[j,i]-q[j-1,i]*rho[j-1,i]<0
                fract1=1-q[j,i]*rho[j,i]/(q[j-1,i]*rho[j-1,i])
                push!(fract,fract1 )
                xm=0.5*(x[j,i]+x[j,i])
                ym=0.5*(y[j,i]+y[j,i])
                pm=0.5*(p[j,i]+p[j,i])
                push!(ix,xm)
                push!(iy,ym)
                push!(ip,pm)
                push!(itime,j)
                q1=q[j-1,i]
                for k=j-1:-1:2
                    if q[k,i]*rho[k,i]-q[k-1,i]*rho[k-1,i]>0
                        q1=q1*q[k-1,i]*rho[k-1,i]/(q[k,i]*rho[k,i])
                    end
                end
                push!(qs,q1)
                old=old*(1-fract1)
                for k=1:size(fract)[1]-1
                    fract[k]=fract[k]*(1-fract1)
                end
            end
        end
        for j=nt:-1:2
            if q[j,i]*rho[j,i]-q[j-1,i]*rho[j-1,i]>0
                qold=qold*q[j-1,i]*rho[j-1,i]/(q[j,i]*rho[j,i])
            end
        end
        xm=0.5*(x[nt,i]+x[nt,i])
        ym=0.5*(y[nt,i]+y[nt,i])
        pm=0.5*(p[nt,i]+p[nt,i])
        inty=Int(trunc(ym)+1)
        intx=Int(trunc(xm+180.0)+1)
        intp=Int(trunc(pm/dp)+1)
        oldM+=old*q[1,i]*rho[1,i]
        evapM+=sum(fract)*q[1,i]*rho[1,i]
        if intx>=1 && intx<=360 && inty>=1 && inty<=90 &&
            intp>=1 && intp<=20
            #old2d[intx,inty,intp]=old2d[intx,inty,intp]+old*q[1,i]
            #cold2d[intx,inty,intp]=cold2d[intx,inty,intp]+size(fract)[1]
        end
        #println(size(fract)[1])
        #println(ix)
        #println(iy)
        #println(ip)
        for k=1:size(fract)[1]
            xm=ix[k]
            ym=iy[k]
            pm=ip[k]
            inty=Int(trunc(ym)+1)
            intx=Int(trunc(xm+180.0)+1)
            intp=Int(trunc(pm/dp)+1)
            #println(intx," ",inty," ",intp)
            if intx>=1 && intx<=360 && inty>=1 && inty<=90 &&
                intp>=1 && intp<=20
                if itime[k]==nt
                    src[intx,inty,intp,k]=src[intx,inty,intp,k]+fract[k]*q[1,i]
                    csrc[intx,inty,intp,k]=csrc[intx,inty,intp,k]+1
                    old2d[intx,inty,intp]=old2d[intx,inty,intp]+old*q[1,i]
                    cold2d[intx,inty,intp]=cold2d[intx,inty,intp]+1
                else
                    src[intx,inty,intp,k]=src[intx,inty,intp,k]+fract[k]*q[1,i]
                    csrc[intx,inty,intp,k]=csrc[intx,inty,intp,k]+1
                end
            end
        end
        #println(old+sum(fract), " ",old)
        
        #println("fract=",fract)
        #println("ix=",ix)
        #println(sum(fract)+old)
        #exit(0)
    end
    #exit(0)
    return src, csrc, old2d, cold2d, oldM, evapM
end
