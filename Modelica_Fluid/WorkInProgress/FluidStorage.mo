package FluidStorage 
  
connector FluidPort_ArrayIcon 
    "Fluid connector with icon suited for an array of FluidPorts" 
    import Modelica_Fluid;
  extends Modelica_Fluid.Interfaces.FluidPort;
  annotation (defaultComponentName="ports",
              Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
             fillColor=69)), Rectangle(extent=[-100, 100; 100, -100], style(color=16,
             fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
       Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
            fillColor=69)), Rectangle(extent=[-100, 100; 100, -100], style(color=16,
            fillColor=69))));
end FluidPort_ArrayIcon;
  
model TankAttachment "Equations to attach pipe at tank" 
      replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
    
    Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
    annotation (extent=[-10,-112; 10,-92],    rotation=90);
   // Real mXi_flow;
    parameter Boolean onlyInFlow = false 
      "= true, if flow only into the tank (e.g. top ports)" 
                                                          annotation(Evaluate=true);
    parameter SI.Diameter pipeDiameter=0.0 "Inner (hydraulic) pipe diameter" 
                                                                         annotation(Dialog(enable=not onlyInFlow));
    parameter SI.Height pipeHeight "Height of pipe";
    parameter Medium.SpecificEnthalpy h_start 
      "Start value of specific enthalpy (used as enthalpy for topPorts if back flow)";
    parameter Medium.MassFraction X_start[Medium.nX] 
      "Start value of mass fractions m_i/m";
    parameter SI.Height level_start "Start value of tank level";
    
    parameter SI.AbsolutePressure p_ambient "Tank surface pressure" 
                                                                annotation(Dialog);
    input SI.Height level "Actual tank level" annotation(Dialog);
    input Medium.SpecificEnthalpy h "Actual specific enthalpy of fluid in tank"
                                                                annotation(Dialog);
    input Medium.Density d "Actual specific density of fluid in tank" 
                                                      annotation(Dialog);
    input Medium.MassFraction Xi[Medium.nXi] 
      "Actual mass fractions of fluid in tank"                  annotation(Dialog);
    parameter Real k_small(min=0) = 1e-5 
      "Small regularization range if tank level is below bottom_height or side_height; k_small = 0 gives ideal switch"
              annotation(Evaluate=true);
    parameter Real s_start = 0;
    
    output Medium.EnthalpyFlowRate H_flow 
      "= port.H_flow (used to transform vector of connectors in vector of Real numbers)";
    output Medium.MassFlowRate m_flow 
      "= port.m_flow (used to transform vector of connectors in vector of Real numbers)";
    output Medium.MassFlowRate mXi_flow[Medium.nXi] 
      "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
    
  annotation (Documentation(info="<html>
<p>
This component contains the equations that attach the pipe
to the tank. The main reason to introduce this component is
that Dymola has currently limitations for connector arrays
when the dimension is zero. Without this utility component
it would not be possible to set, e.g., n_topPorts to zero.
</p>
 
 
</html>"),
         Icon(Rectangle(extent=[-100,0; 100,-100], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255})), Text(
        extent=[-122,48; 132,6],
        style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1),
        string="%name")));
  Modelica_Fluid.Ambient ambient;
  protected 
  parameter SI.Area pipeArea = Modelica.Constants.pi*(pipeDiameter/2)^2;
  parameter Medium.MassFlowRate m_flow_nominal = 1 
      "Nominal mass flow rate used for scaling (has only an effect if k_small > 0)";
  parameter Medium.AbsolutePressure p_nominal = p_ambient;
  SI.Length aboveLevel = level - pipeHeight;
  Boolean m_flow_out(start=true,fixed=true) "true= massflow out of tank";
  Real s(start=s_start) 
      "path parameter of parameterized curve description (either m_flow/m_flow_nominal or (port.p-p_ambient)/p_ambient)";
equation 
  H_flow = port.H_flow;
  m_flow = port.m_flow;
  mXi_flow = port.mXi_flow;
    
  if onlyInFlow then
     m_flow_out = false "Dummy value in this if branch";
     port.p = p_ambient;
     /* flow should never out of the port. However, if this occurs in a 
        small time interval (e.g. during initialization), the start values of
        h and X are provided, since otherwise there is a singular
        system 
     */
     port.H_flow = semiLinear(port.m_flow, port.h, h_start);
     port.mXi_flow = semiLinear(port.m_flow, port.Xi, X_start[1:Medium.nXi]);
     assert(port.m_flow > -1e-6, "Mass flows out of tank via topPort. This indicates a wrong model");
     s = 0;
  else
     port.H_flow = semiLinear(port.m_flow, port.h, h);
     port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
      
/* Original equations from Poschlad/Remelhe:
*/
     s = 0;
     m_flow_out = (pre(m_flow_out) and not port.p>p_ambient) or (port.m_flow < -1e-6);
      
     if (aboveLevel > 0) then
       port.p = aboveLevel*ambient.g*d + p_ambient - smooth(2,noEvent(if m_flow < 0 then m_flow^2/(2*d*pipeArea^2) else 0));
     else
       if pre(m_flow_out) then
          m_flow = 0;
       else
          port.p = p_ambient;
       end if;
     end if;
      
/* Martin Otter: The following equations are a declarative form 
   (parameterized curve description) of the above equations and
   should theoretically work better. However, some examples with
   IF97 water fail, whereas the above works. Therefore, not used. 
       Add the following text to OpenTank, once the initialization
   with this solution works for Modelica_Fluid.Examples.Tanks.ThreeOpenTanks:
 
OpenTank:
<p>
The situation when the tank level is below bottom_heights[i] or side_heights[i]
is handeled properly. Details are described
<a href="Modelica:Modelica_Fluid.Utilities.TankAttachment">here</a>
</p> 
 
TankAttachment:
<p>
If a bottom or side connector is above the actual tank level, the
following characteristic is used to compute the mass flow rate port.m_flow
from the connector to the tank and the absolute pressure port.p
in the port:
</p>
 
<img src="../Images/Components/Tank_PipeAboveTankLevel.png">   
 
 
     m_flow_out = s <= 0;
     if aboveLevel >= 0 then
        m_flow = m_flow_nominal*s "equation to compute s, which is a dummy in this branch";
        port.p - p_ambient = aboveLevel*fluidOptions.g*d  -
                             smooth(2,if m_flow_out then s*abs(s)*m_flow_nominal^2/(2*d*pipeArea^2) else k_small*m_flow_nominal*s);
     else
        m_flow = (if m_flow_out then k_small*p_nominal else m_flow_nominal)*s;
        port.p - p_ambient = (if m_flow_out then p_nominal else k_small*m_flow_nominal)*s;
     end if;
*/
      
  end if;
    
  /*
  More precise equations (introduce them later; need to transform
  them from energy balance form 1 to form 2):
 
  Momentum balance:
  (integrated momentum equation for frictionless fluid with density that is
   independent of the level, i.e., the unsteady Bernoulli equation for incompressible fluid)
  v_level = der(level);
  v = -port.m_flow/(rho*A_outlet);
  level*der(v_level) + (v^2 - v_level^2)/2 - g*level + (p - p_ambient)/rho = 0; or
  rho*level*der(v_level) + rho*(v^2 - v_level^2)/2 - rho*g*level + (p - p_ambient) = 0;
 
  Energy balance:
  Potential energy: E_pot = integ(dm*g*s)
                          = g*integ(rho*A*s*ds)
                          = g*rho*A*z^2/2
  Kinetic energy  : E_kin = integ(dm*v^2/2)
                          = integ(rho*A*v^2/2*ds)
                          = rho*A*v^2/2*integ(ds)
                          = rho*A*v^2/2*z
                          = M*v^2/2
  E = U + M*g*z/2 + M*v_level^2/2
  der(E) = port.H_flow + port.m_flow*v^2/2 - p_ambient*area*der(level)
*/
    
end TankAttachment;
end FluidStorage;
