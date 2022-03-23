
begin

infiledir   = "./data/"
infilename  = "annmean-pi0400co21f-1501-1700-3dlvl.nc"

newfiledir  = "./thermocline/"
newfilename = "max-curvature-fin-thermocline-pi0400co21f.nc"

;==========================================
print("Reading file...")

data    = addfile(infiledir + infilename, "r")

temp  = data->templvl	;(1,70,385,360)
depth = data->depth	;(70) 
lon   = data->plon	;(385,360)
lat   = data->plat

;===========================================
print("Defining new valuables...")

thermocline = new((/385,360/), double, 9.96920996838687e+36)
thermocline!0 = "y"
thermocline!1 = "x"
thermocline@long_name = "Thermocline Depth"
thermocline@units = "m"
thermocline@_FillValue = 9.96920996838687e+36
thermocline@cell_measures = "area: parea"

;===========================================
print("Calculating thermocline depth...")

h = 0

do ny = 0,384
    do nx = 0,359

	;--- h 作为储存该格点最大的温度梯度
        h = 0
        do nz = 1,27
	   
            ;--- 跳过陆地区域格点
            if(ismissing(temp(0,nz,ny,nx)))then
                break
            end if

 	    ;--- 计算温度梯度
            k = (temp(0,nz-1,ny,nx)-temp(0,nz,ny,nx)) / (depth(nz)-depth(nz-1))
            
	    if(k.gt.h)then
                thermocline(ny,nx) =(/depth(nz-1)/)
                h := k
            end if
        end do

    end do
end do

do ny = 0,384
    do nx = 0,359

	h = 0	

	;--- 200 ~ 500m 的深水层
	do nz = 19,27

	    if(ismissing(temp(0,nz,ny,nx)))then
	        break
	    end if

	    k = (temp(0,nz-1,ny,nx)-temp(0,nz,ny,nx)) / (depth(nz)-depth(nz-1))
	    
	    ;--- 条件筛选，若不筛选则把 if 和 end if 注释掉
	    if(k.ge.0.05)then

		;--- 比较最大的温度梯度，并储存其所在的深度    
	        if(k.gt.h)then
	  	    thermocline(ny,nx) =(/depth(nz-1)/)
	            h := k
		end if 

	    end if
	end do
	
	;--- 200m 以上的浅水层	
	do nz = 1,18	    
            if(ismissing(temp(0,nz,ny,nx)))then
                break
            end if
            k = (temp(0,nz-1,ny,nx)-temp(0,nz,ny,nx)) / (depth(nz)-depth(nz-1))
	    
	    if(k.ge.0.2)then

                if(k.gt.h)then
                    thermocline(ny,nx) =(/depth(nz-1)/)
                    h := k
                end if
	    
	    end if
	end do


    end do
end do

;============================================
print("Creating new file...")

system("rm "+ newfiledir + newfilename)
newfile = addfile(newfiledir + newfilename,"c")

;============================================
print("Writing new file...")

fileAtt               = True
fileAtt@tittle        = "Thermocline"
fileAtt@Conventions   = "None"
fileAtt@creation_date = systemfunc("date")
fileattdef(newfile,fileAtt)

dimNames = (/"time","depth","y","x"/)
dimSizes = (/1,70,385,360/)
dimUnlim = (/False,False,False,False/)

filedimdef(newfile,dimNames,dimSizes,dimUnlim)

filevardef(newfile,"thermocline",typeof(thermocline),(/"y","x"/))
filevarattdef(newfile,"thermocline",thermocline)
newfile->thermocline = (/thermocline/)

newfile->templvl = temp
newfile->plon = lon
newfile->plat = lat
newfile->depth = depth

;============================================
print("Program Done!")
end
