
begin

;=============== Name of Files
filename1 = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"    
filename2 = "surfdata_1.9x2.5_simyr1850_c091108.nc"

newfile1  = "new_USGS-gtopo30_1.9x2.5_remap_c050602.nc"
newfile2  = "new_surfdata_1.9x2.5_simyr1850_c091108.nc"

;=============== Create Newfile
print("Creating new file...")

system("cp "+ filename1 + " ./"+ newfile1)
system("cp "+ filename2 + " ./"+ newfile2)

;=============== Read File
print("Reading file...")

data1 = addfile(newfile1,"w")
data2 = addfile(newfile2,"w")

PHIS        = data1->PHIS
PCT_GLACIER = data2->PCT_GLACIER
PCT_PFT     = data2->PCT_PFT

;======= Modify PHIS 
print("Calculating PHIS...")

do nj = 82, 90 
    do ni = 125, 133
	if(PHIS(nj,ni).lt.15000)then        
		PHIS(nj,ni) = 0
	end if
	if(PHIS(nj,ni).ge.15000)then
		PHIS(nj,ni) = PHIS(nj,ni) / 2
    	end if
    end do
end do

data1->PHIS = PHIS

;======= Modify Surface Data
print("Calculating Surface Data...")

do nj = 77, 94 
    do ni = 100, 140
        if(PCT_GLACIER(nj,ni).ne.100)then
		PCT_PFT(11,nj,ni)  = PCT_PFT(11,nj,ni) + PCT_GLACIER(nj,ni) / 2
		PCT_PFT(0,nj,ni)   = PCT_PFT(0,nj,ni)  + PCT_GLACIER(nj,ni) / 2
                PCT_GLACIER(nj,ni) = 0
        end if
    end do
end do

do nj = 77, 94  
    do ni = 100, 140
        if((PCT_GLACIER(nj+1,ni).eq.0).or.(PCT_GLACIER(nj-1,ni).eq.0))then    
	        PCT_PFT(11,nj,ni)  = PCT_PFT(11,nj,ni) + PCT_GLACIER(nj,ni) / 4
		PCT_PFT(0,nj,ni)   = PCT_PFT(0,nj,ni)  + PCT_GLACIER(nj,ni) / 4
	 	PCT_GLACIER(nj,ni) = PCT_GLACIER(nj,ni) / 2
        end if

	if((PCT_GLACIER(nj,ni+1).eq.0).or.(PCT_GLACIER(nj,ni-1).eq.0))then 
                PCT_PFT(11,nj,ni)  = PCT_PFT(11,nj,ni) + PCT_GLACIER(nj,ni) / 4
		PCT_PFT(0,nj,ni)   = PCT_PFT(0,nj,ni)  + PCT_GLACIER(nj,ni) / 4
		PCT_GLACIER(nj,ni) = PCT_GLACIER(nj,ni) / 2
        end if
    end do
end do

data2->PCT_GLACIER = PCT_GLACIER
data2->PCT_PFT     = PCT_PFT

;================
print("Promgram Done!")

end

