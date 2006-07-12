package ControlVolumes 
model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
    import Modelica_Fluid.Types;
  extends Modelica_Fluid.Interfaces.Records.PartialInitializationParameters;
    
  replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    
  parameter SI.Volume V "Volume";
    
  Modelica_Fluid.Interfaces.Ports.FluidPort_a port(
    redeclare package Medium = Medium) "Fluid port" 
    annotation (extent=[-10, -10; 10, 10], rotation=0);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
    annotation (extent=[-10,90; 10,110]);
    
  Medium.BaseProperties medium(preferredMediumStates=true,
              p(start=p_start), T(start=T_start),
              h(start=h_start), Xi(start= X_start[1:Medium.nXi]));
  SI.Energy U "Internal energy of fluid";
  SI.Mass m "Mass of fluid";
  SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
  annotation (
   Icon(
      Ellipse(extent=[-100, 100; 100, -100], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=3,
          fillColor=68,
          rgbfillColor={170,213,255})),
      Text(extent=[-150,-100; 150,-150], string="%name")),
                         Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (= medium.h) is always identical to the specific enthalpy
in the port (port.h = medium.h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport. This splitting is only possible under certain assumptions.
</p>
</html>"),
    Diagram);
equation 
  // medium properties set by port values
    port.p = medium.p;
    port.h = medium.h;
    port.Xi = medium.Xi;
    thermalPort.T = medium.T;
    
  // Total quantities
     m    = V*medium.d "Total Mass";
     mXi = m*medium.Xi "Independent component masses";
     U    = m*medium.u "Internal energy";
    
  // Mass and energy balance
     der(m)    = port.m_flow "Total mass balance";
     der(mXi)  = port.mXi_flow "Independent component mass balance";
     der(U)    = port.H_flow + thermalPort.Q_flow "Energy balance";
    
initial equation 
  // Initial conditions
  if initType == Types.Init.NoInit then
    // no initial equations
  elseif initType == Types.Init.InitialValues then
    if not Medium.singleState then
       medium.p = p_start;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  elseif initType == Types.Init.SteadyState then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    der(medium.h) = 0;
    der(medium.Xi) = zeros(Medium.nXi);
  elseif initType == Types.Init.SteadyStateHydraulic then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  else
    assert(false, "Unsupported initialization option");
  end if;
end PortVolume;
end ControlVolumes;
