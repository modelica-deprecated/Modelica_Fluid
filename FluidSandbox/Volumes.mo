within FluidSandbox;
package Volumes "Generic mixing volume and other volume type components" 
  extends Icons.VariantLibrary;
  
  model Volume "Volume with inlet and outlet ports (flow reversal is allowed)" 
    
    import Modelica_Fluid.Types;
    
    extends Interfaces.PartialComponent;
    
    extends FluidInterface.PartialLumpedVolume(medium(
      preferredMediumStates=true,
      p(start=p_start),
      h(start=h_start),
      T(start=T_start),
      Xi(start=X_start[1:Medium.nXi])));
    
    // Constant volume
    parameter SI.Volume V "Volume";
    
  // Extensive properties
    SI.Energy U "Internal energy of fluid";
    SI.Mass m "Mass of fluid";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    
  //Initialization  
    parameter Types.Init.Temp initType=Types.Init.NoInit 
      "Initialization option" 
    annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start=Medium.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start=true "= true, use T_start, otherwise h_start"
    annotation(Dialog(tab = "Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
              Medium.temperature_phX(
            p_start,
            h_start,
            X_start) "Start value of temperature" 
    annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
        Medium.specificEnthalpy_pTX(
            p_start,
            T_start,
            X_start) else Medium.h_default "Start value of specific enthalpy" 
    annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    
  equation 
  // Extensive quantities
    m = V*medium.d;
    mXi = m*medium.Xi;
    U = m*medium.u;
    
  // Mass and energy balance
    der(m) = m_flow_net;
    der(mXi) = mXi_flow_net;
    der(U) = H_flow_net + heatPort.Q_flow;
    
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
    annotation (
      Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-144,178; 146,116], string="%name"),
        Text(
          extent=[-130,-108; 144,-150],
          style(color=0),
          string="V=%V")),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. The thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram);
  end Volume;
end Volumes;
