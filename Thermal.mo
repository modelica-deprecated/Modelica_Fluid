package Thermal "Components to model the thermal behavior of pipe walls" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  
model WallConstProps 
    "Pipe wall with capacitance, assuming 1D heat conduction and constant material properties" 
  constant Real pi=Modelica.Constants.pi;
  parameter Integer n(min=1)=1 "Segmentation perpendicular to heat conduction";
  parameter SI.Diameter a_inner "Inner cross section area";
  parameter SI.Length a_outer "Outer cross section area";
  parameter SI.Length length "Pipe length";
  parameter SI.Length area_h "Heat transfer area";
  parameter SI.Length s=sqrt(a_outer/pi)-sqrt(a_inner/pi) 
      "Wall thickness, default expression for circular pipe";
  SI.Temperature[n] Tb(each start=T_start+0.5*dT);
  SI.Temperature[n] Ta(each start=T_start-0.5*dT);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_a 
      "Thermal port" 
    annotation (extent=[-20,40; 20,60]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_b 
      "Thermal port" 
    annotation (extent=[-20,-40; 20,-60]);
  parameter SI.Density d_wall "Density of wall material";
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material";
  parameter SI.ThermalConductivity k_wall 
      "Thermal conductivity of wall material";
  parameter SI.Temperature T_start "Start value for wall temperature";
  parameter SI.Mass[n] m=ones(n)*(a_outer-a_inner)*length*d_wall/n 
      "Distribution of wall mass";
  parameter Types.Init.Temp initType=Types.
       Init.NoInit "Initialization option" 
   annotation(Evaluate=true, Dialog(tab = "Initialization"));
  parameter Boolean initWall_steadyState=false 
      " = true, wall is initialized in steady state";
  SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer) 
      "Wall temperature";
  parameter SI.Temperature dT "Start value for port_b.T - port_a.T";
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
   thermalPort_a[i].Q_flow=k_wall/s*(Ta[i]-T[i])*area_h/n;
   thermalPort_b[i].Q_flow=k_wall/s*(Tb[i]-T[i])*area_h/n;
  end for;
  Ta=thermalPort_a.T;
  Tb=thermalPort_b.T;
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
Simple model of circular (or any other closed shape) wall to be used for pipe (or duct) models. Heat conduction is regarded one dimensional, capacitance is lumped at the arithmetic mean temperature. The spatial discretization (parameter <tt>n</tt>) is meant to correspond to a connected fluid model discretization.
</html>"));
end WallConstProps;
  
  annotation (Documentation(info="<html>
 
</html>"));
end Thermal;
