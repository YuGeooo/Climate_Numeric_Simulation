
begin

;=============== Name of Files

filename = "surfdata_1.9x2.5_simyr1850_c091108.nc"    
newfile  = "new_surfdata_1.9x2.5_simyr1850_c091108.nc"

;=============== Create Newfile
print("Creating new file...")

system("cp "+ filename + " ./"+ newfile)

;=============== Read File
print("Reading file...")

data = addfile(newfile,"w")
PCT_PFT = data->PCT_PFT

;======= Transform SaHara Desert into Evergreen Tropical Tree 
print("Calculating...")

do nj = 56, 64 
    do ni = 0, 24
        PCT_PFT(4,nj,ni) = PCT_PFT(0,nj,ni) + PCT_PFT(4,nj,ni)
        PCT_PFT(0,nj,ni) = 0
    end do
end do

data->PCT_PFT = PCT_PFT

;================
print("Promgram Done!")

end

