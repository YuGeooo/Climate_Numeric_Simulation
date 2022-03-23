
begin

infiledir   = "./data/"
infilename  = "annmean-pi0400co21f-1501-1700-3dlvl.nc"

newfiledir  = "./thermocline/"
newfilename = "new-new-thermocline-pi0400co21f.nc"

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

do ny = 0,384
    do nx = 0,359

	do nz = 19,69
	    if(ismissing(temp(0,nz,ny,nx)))then
	        break
	    end if
	    k = (temp(0,nz-1,ny,nx)-temp(0,nz,ny,nx)) / (depth(nz)-depth(nz-1))
	    if(k.gt.0.05)then
		thermocline(ny,nx) =(/depth(nz-1)/)
	        break
	    end if
	end do

	do nz = 1,18
            if(ismissing(temp(0,nz,ny,nx)))then
                break
            end if
	    k = (temp(0,nz-1,ny,nx)-temp(0,nz,ny,nx)) / (depth(nz)-depth(nz-1))
	    if(k.gt.0.2)then
		thermocline(ny,nx) = (/depth(nz-1)/)
		break
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
