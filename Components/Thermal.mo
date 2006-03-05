package Thermal 
model WallConstProps 
    "Pipe wall with capacitance, assuming ideal 1D-conduction and constant material properties" 
  parameter Integer n(min=1)=1 "Pipe segmentation";
  parameter SI.Diameter a_inner "Inner cross section area";
  parameter SI.Length a_outer "Outer cross section area";
  parameter SI.Length length "Pipe length";
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_a 
      "Thermal port" 
    annotation (extent=[-20,40; 20,60]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_b 
      "Thermal port" 
    annotation (extent=[-20,-40; 20,-60]);
  parameter SI.Density d_wall "Density of wall material";
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material";
  parameter SI.Temperature T_start "Start value for wall temperature";
  parameter SI.Mass[n] m=ones(n)*(a_outer-a_inner)*length*d_wall/n "Wall mass";
  parameter Boolean initWall_steadyState=false 
      " = true, wall is initialized in steady state";
  SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer) 
      "Wall temperature";
initial equation 
  if initWall_steadyState then
    der(T)=zeros(n);
  else
   T=ones(n)*T_start;
  end if;
equation 
    
  for i in 1:n loop
   assert(m[i]>0, "Wall has negative dimensions");
   c_wall*m[i]*der(T[i]) = thermalPort_a[i].Q_flow + thermalPort_b[i].Q_flow;
  end for;
  //assuming ideal heat conduction perpendicular to fluid flow, conduction in remaining two dimensions is neglected
  thermalPort_a.T=T;
  thermalPort_b.T=T;
    annotation (Icon(      Rectangle(extent=[-100,40; 100,-40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)), Text(
        extent=[-82,18; 76,-18],
        string="%name",
        style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=7))), Documentation(revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>", info="<html>
Simple model of circular (or any other closed shape) wall to be used for pipe (or duct) models. Ideal radial heat conductance is assumed, conductance in other directions is neglected.
</html>"));
end WallConstProps;
  
end Thermal;
