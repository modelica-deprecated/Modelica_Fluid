model TestPressureDropDerivatives 
  "Test that PressureDrop components can be differentiated" 
  import Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent;
  extends Modelica.Icons.Example;
  parameter SI.Diameter D_a = 0.1 "Inner diameter of pipe at port_a";
  parameter SI.Diameter D_b = 0.2 "Inner diameter of pipe at port_b";
  parameter QuadraticTurbulent.LossFactorData data=
                  QuadraticTurbulent.LossFactorData.suddenExpansion(D_a, D_b) 
    "Loss factors for both flow directions";
  parameter Real dp_small = 0.1;
  
  Real d_a;
  Real d_b;
  Real eta_a;
  Real eta_b;
  Real dp;
  Real m_flow1;
  Real m_flow2;
  Real der_m_flow1;
  // Real der_m_flow2;
  
  annotation (experiment(StopTime=3), experimentSetupOutput,
    Documentation(info="<html>
 
</html>"));
equation 
  dp = time - 1;
  d_a = 0.1 + time;
  d_b = 0.2 + 0.5*time;
  eta_a = 0.1*d_a;
  eta_b = 0.4*d_b;
  
  m_flow1 = QuadraticTurbulent.massFlowRate_dp(dp, d_a, d_b, data, dp_small);
  m_flow2 = QuadraticTurbulent.massFlowRate_dp_and_Re(dp, d_a, d_b, eta_a, eta_b, data);
  
  der_m_flow1 = der(m_flow1);
  // der_m_flow2 = der(m_flow2);
end TestPressureDropDerivatives;
