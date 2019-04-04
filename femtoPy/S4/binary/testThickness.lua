function printFluxOrders(file,layer,offset)
   G=S:GetGList()
   P=S:GetPowerFluxByOrder(layer,offset)
   for i=1,numG,1
   do
      file:write(G[i][1]," ",G[i][2]," ")
      for j=1,4,1
      do
	 file:write(P[i][j]," ")
      end
      file:write("\n")
   end
   file:write("\n")
end

input=io.open("init.dat","r")
pcall(loadstring(input:read("*all")))
trans=io.open(trans,"w")
refl=io.open(refl,"w")

for t=tmin,tmax,tstep
do
   S=S4.NewSimulation()
   
   S:SetLattice({1,0},{0,0})
   S:SetNumG(numG)
   
   S:AddMaterial("Vacuum",{1,0})
   --S:AddMaterial("Gmat",{ngr,ngi})
   S:AddMaterial("Gmat",{13,0})
   --S:AddMaterial("Substrate",{nsr,nsi})
   S:AddMaterial("Substrate",{1,0})
   
   S:AddLayer("Vac",0,"Vacuum")
   --S:AddLayer("Grat",t,"Vacuum")
   S:AddLayer("Grat",t,"Vacuum")
   S:SetLayerPatternRectangle("Grat","Gmat",{center,0},0,{halfWidth,0})
   S:AddLayer("Sub",0,"Substrate")

   S:SetExcitationPlanewave({theta,0},{s_amp,s_phase},{p_amp,p_phase})
   S:SetFrequency(f)

   printFluxOrders(trans,"Sub",0)
   printFluxOrders(refl,"Vac",0)
end
   trans.close()
   refl.close()
