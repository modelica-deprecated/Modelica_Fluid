within Modelica_Fluid;
package Thermal "Components to model the thermal behavior of pipe walls"
  extends Modelica_Fluid.Icons.VariantLibrary;

model WallConstProps
    "Pipe wall with capacitance, assuming 1D heat conduction and constant material properties"
  parameter Integer n(min=1)=1 "Segmentation perpendicular to heat conduction";
//Geometry
  parameter SI.Length s "Wall thickness";
  parameter SI.Area area_h "Heat transfer area";
//Material properties
  parameter SI.Density d_wall "Density of wall material";
  parameter SI.SpecificHeatCapacity c_wall
      "Specific heat capacity of wall material";
  parameter SI.ThermalConductivity k_wall
      "Thermal conductivity of wall material";
  parameter SI.Mass[n] m=fill(d_wall*area_h*s/n,n) "Distribution of wall mass";
//Initialization
  parameter Types.Init initType=Types.
        Init.NoInit "Initialization option" 
    annotation(Evaluate=true);
  parameter SI.Temperature T_start "Wall temperature start value";
  parameter SI.Temperature dT "Start value for port_b.T - port_a.T";
//Temperatures
  SI.Temperature[n] Tb(each start=T_start+0.5*dT);
  SI.Temperature[n] Ta(each start=T_start-0.5*dT);
  SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer)
      "Wall temperature";
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_a
      "Thermal port" 
    annotation (Placement(transformation(extent={{-20,40},{20,60}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_b
      "Thermal port" 
    annotation (Placement(transformation(extent={{-20,-40},{20,-60}}, rotation=
              0)));

initial equation
  if initType == Types.Init.SteadyState or initType == Types.Init.SteadyStateHydraulic then
    der(T) = zeros(n);
  elseif initType == Types.Init.InitialValues then
    T = ones(n)*T_start;
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
    annotation (Icon(graphics={Rectangle(
            extent={{-100,40},{100,-40}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Forward), Text(
            extent={{-82,18},{76,-18}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Forward,
            textString=
               "%name")}),  Documentation(revisions="<html>
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
