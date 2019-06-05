/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2014 Universidad de la República - Uruguay
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Matias Richart <mrichart@fing.edu.uy>
 */

/**
 * This example program is designed to illustrate the behavior of three
 * power/rate-adaptive WiFi rate controls; namely, ns3::ParfWifiManager,
 * ns3::AparfWifiManager and ns3::RrpaaWifiManager.
 *
 * The output of this is typically two plot files, named throughput-parf.plt
 * (or throughput-aparf.plt, if Aparf is used) and power-parf.plt. If
 * Gnuplot program is available, one can use it to convert the plt file
 * into an eps file, by running:
 * \code{.sh}
 *   gnuplot throughput-parf.plt
 * \endcode
 * Also, to enable logging of rate and power changes to the terminal, set this
 * environment variable:
 * \code{.sh}
 *   export NS_LOG=PowerAdaptationDistance=level_info
 * \endcode
 *
 * This simulation consist of 2 nodes, one AP and one STA.
 * The AP generates UDP traffic with a CBR of 54 Mbps to the STA.
 * The AP can use any power and rate control mechanism and the STA uses
 * only Minstrel rate control.
 * The STA can be configured to move away from (or towards to) the AP.
 * By default, the AP is at coordinate (0,0,0) and the STA starts at
 * coordinate (5,0,0) (meters) and moves away on the x axis by 1 meter every
 * second.
 *
 * The output consists of:
 * - A plot of average throughput vs. distance.
 * - A plot of average transmit power vs. distance.
 * - (if logging is enabled) the changes of power and rate to standard output.
 *
 * The Average Transmit Power is defined as an average of the power
 * consumed per measurement interval, expressed in milliwatts.  The
 * power level for each frame transmission is reported by the simulator,
 * and the energy consumed is obtained by multiplying the power by the
 * frame duration.  At every 'stepTime' (defaulting to 1 second), the
 * total energy for the collection period is divided by the step time
 * and converted from dbm to milliwatt units, and this average is
 * plotted against time.
 *
 * When neither Parf, Aparf or Rrpaa is selected as the rate control, the
 * generation of the plot of average transmit power vs distance is suppressed
 * since the other Wifi rate controls do not support the necessary callbacks
 * for computing the average power.
 *
 * To display all the possible arguments and their defaults:
 * \code{.sh}
 *   ./waf --run "power-adaptation-distance --help"
 * \endcode
 *
 * Example usage (selecting Aparf rather than Parf):
 * \code{.sh}
 *   ./waf --run "power-adaptation-distance --manager=ns3::AparfWifiManager --outputFileName=aparf"
 * \endcode
 *
 * Another example (moving towards the AP):
 * \code{.sh}
 *   ./waf --run "power-adaptation-distance --manager=ns3::AparfWifiManager --outputFileName=aparf --stepsSize=-1 --STA1_x=200"
 * \endcode
 *
 * To enable the log of rate and power changes:
 * \code{.sh}
 *   export NS_LOG=PowerAdaptationDistance=level_info
 * \endcode
 */

#include "ns3/gnuplot.h"
#include "ns3/command-line.h"
#include "ns3/config.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "ns3/log.h"
#include "ns3/yans-wifi-helper.h"
#include "ns3/ssid.h"
#include "ns3/mobility-helper.h"
#include "ns3/internet-stack-helper.h"
#include "ns3/ipv4-address-helper.h"
#include "ns3/packet-sink-helper.h"
#include "ns3/on-off-helper.h"
#include "ns3/yans-wifi-channel.h"
#include "ns3/wifi-net-device.h"
#include "ns3/wifi-mac.h"
#include "ns3/wifi-mac-header.h"
#include "ns3/mobility-model.h"
#include "ns3/propagation-loss-model.h"

using namespace ns3;
using namespace std;

NodeContainer wifiApNodes;
NodeContainer wifiStaNodes;

NetDeviceContainer wifiApDevices;
NetDeviceContainer wifiStaDevices;
NetDeviceContainer wifiDevices;

MobilityHelper mobility;

double maxPower = 17;
double minPower = 0;

double frequency = 5.150e9;
double C = 299792458.0;
double m_lambda = C / frequency;
double m_systemLoss = 1.0;
double m_minLoss = 0;
double m_pi = 3.141592;

double max(double a, double b){
  return a>b?a:b;
}

struct point 
{
    double x,y;
};

double norm (point p) // get the norm of a vector
{
    return pow(pow(p.x,2)+pow(p.y,2),.5);
}

point trilateration(point point1, point point2, point point3, double r1, double r2, double r3) {
    point resultPose;
    //unit vector in a direction from point1 to point 2
    double p2p1Distance = pow(pow(point2.x-point1.x,2) + pow(point2.y-   point1.y,2),0.5);
    point ex = {(point2.x-point1.x)/p2p1Distance, (point2.y-point1.y)/p2p1Distance};
    point aux = {point3.x-point1.x,point3.y-point1.y};
    //signed magnitude of the x component
    double i = ex.x * aux.x + ex.y * aux.y;
    //the unit vector in the y direction. 
    point aux2 = { point3.x-point1.x-i*ex.x, point3.y-point1.y-i*ex.y};
    point ey = { aux2.x / norm (aux2), aux2.y / norm (aux2) };
    //the signed magnitude of the y component
    double j = ey.x * aux.x + ey.y * aux.y;
    //coordinates
    double x = (pow(r1,2) - pow(r2,2) + pow(p2p1Distance,2))/ (2 * p2p1Distance);
    double y = (pow(r1,2) - pow(r3,2) + pow(i,2) + pow(j,2))/(2*j) - i*x/j;
    //result coordinates
    double finalX = point1.x+ x*ex.x + y*ey.x;
    double finalY = point1.y+ x*ex.y + y*ey.y;
    resultPose.x = finalX;
    resultPose.y = finalY;
    return resultPose;
}

double DoCalcRxPower (double txPowerDbm, Ptr<MobilityModel> a, Ptr<MobilityModel> b){
   /*
    * Friis free space equation:
    * where Pt, Gr, Gr and P are in Watt units
    * L is in meter units.
    *
    *    P     Gt * Gr * (lambda^2)
    *   --- = ---------------------
    *    Pt     (4 * pi * d)^2 * L
    *
    * Gt: tx gain (unit-less)
    * Gr: rx gain (unit-less)
    * Pt: tx power (W)
    * d: distance (m)
    * L: system loss
    * lambda: wavelength (m)
    *
    * Here, we ignore tx and rx gain and the input and output values 
    * are in dB or dBm:
    *
    *                           lambda^2
    * rx = tx +  10 log10 (-------------------)
    *                       (4 * pi * d)^2 * L
    *
    * rx: rx power (dB)
    * tx: tx power (dB)
    * d: distance (m)
    * L: system loss (unit-less)
    * lambda: wavelength (m)
    */
   double distance = a->GetDistanceFrom (b);
   if (distance < 3*m_lambda)
     {
       NS_LOG_UNCOND ("distance not within the far field region => inaccurate propagation loss value");
     }
   if (distance <= 0)
     {
       return txPowerDbm - m_minLoss;
     }
   double numerator = m_lambda * m_lambda;
   double denominator = 16 * m_pi * m_pi * distance * distance * m_systemLoss;
   double lossDb = -10 * log10 (numerator / denominator);
   //NS_LOG_UNCOND ("distance=" << distance<< "m, loss=" << lossDb <<"dB");
   return txPowerDbm - max(lossDb, m_minLoss);
 }

double doCalcDistance(double Pt, double Pr){
  double loss = Pt - Pr;
  double x = pow(10, (loss / -10));
  double xx = (m_lambda * m_lambda) / x;
  double xxx = xx / m_systemLoss / 16 / m_pi / m_pi;
  double distance = sqrt(xxx);
  return distance;
}

void scheduling(int stepsSize, int stepsTime){
  // signal

  // calculate distance
  Ptr<MobilityModel> mobilitySta = wifiStaNodes.Get(0)->GetObject<MobilityModel> ();
  Ptr<MobilityModel> mobilityAp1 = wifiApNodes.Get(0)->GetObject<MobilityModel> ();
  Ptr<MobilityModel> mobilityAp2 = wifiApNodes.Get(1)->GetObject<MobilityModel> ();
  Ptr<MobilityModel> mobilityAp3 = wifiApNodes.Get(2)->GetObject<MobilityModel> ();
  //NS_LOG_UNCOND("sta position: " << mobilitySta->GetPosition());
  Vector ap1_coord = mobilityAp1->GetPosition();
  Vector ap2_coord = mobilityAp2->GetPosition();
  Vector ap3_coord = mobilityAp3->GetPosition();
  //NS_LOG_UNCOND("ap1 position: " << ap1_coord);
  //NS_LOG_UNCOND("ap2 position: " << ap2_coord);
  //NS_LOG_UNCOND("ap3 position: " << ap3_coord);
  double Rx1 = DoCalcRxPower(maxPower, mobilitySta, mobilityAp1);
  double Rx2 = DoCalcRxPower(maxPower, mobilitySta, mobilityAp2);
  double Rx3 = DoCalcRxPower(maxPower, mobilitySta, mobilityAp3);
  //NS_LOG_UNCOND("RX1 POWER: " << Rx1);
  //NS_LOG_UNCOND("RX2 POWER: " << Rx2);
  //NS_LOG_UNCOND("RX3 POWER: " << Rx3);
  double distance1 = doCalcDistance(maxPower, Rx1);
  double distance2 = doCalcDistance(maxPower, Rx2);
  double distance3 = doCalcDistance(maxPower, Rx3);
  //NS_LOG_UNCOND("distance1: " << distance1);
  //NS_LOG_UNCOND("distance2: " << distance2);
  //NS_LOG_UNCOND("distance3: " << distance3);

  point ap1_double_coord = {ap1_coord.x, ap1_coord.y};
  point ap2_double_coord = {ap2_coord.x, ap2_coord.y};
  point ap3_double_coord = {ap3_coord.x, ap3_coord.y};

  point coord = trilateration(ap1_double_coord,
                                          ap2_double_coord,
                                          ap3_double_coord,
                                          distance1,
                                          distance2,
                                          distance3);

  NS_LOG_UNCOND("coord,  x: " << coord.x << " y: " << coord.y);

  // advance position
  Vector sta_coord = mobilitySta->GetPosition();
  sta_coord.x += 1;
  sta_coord.y -= 0.5;
  mobilitySta->SetPosition(sta_coord);

  if(stepsTime < 5) Simulator::Schedule(Seconds(stepsSize), &scheduling, stepsSize, stepsTime + 1);
}

int main (int argc, char *argv[]){
  uint32_t powerLevels = 18;

  uint32_t rtsThreshold = 2346;
  std::string manager = "ns3::ParfWifiManager";
  int ap1_x = 1;
  int ap1_y = 1;
  int ap2_x = 1;
  int ap2_y = 4;
  int ap3_x = 6;
  int ap3_y = 1;
  int sta1_x = 3;
  int sta1_y = 3;
  int packetSize = 5;
  double simuTime = 0.5;
  //uint32_t steps = 200;

  //Define the APs
  wifiApNodes.Create (3);

  //Define the STAs
  wifiStaNodes.Create (1);

  WifiHelper wifi;
  wifi.SetStandard (WIFI_PHY_STANDARD_80211a);
  WifiMacHelper wifiMac;
  YansWifiPhyHelper wifiPhy = YansWifiPhyHelper::Default ();
  YansWifiChannelHelper wifiChannel = YansWifiChannelHelper::Default ();

  wifiPhy.SetChannel (wifiChannel.Create ());

  //Configure the STA node
  wifi.SetRemoteStationManager ("ns3::MinstrelWifiManager", "RtsCtsThreshold", UintegerValue (rtsThreshold));
  wifiPhy.Set ("TxPowerStart", DoubleValue (maxPower));
  wifiPhy.Set ("TxPowerEnd", DoubleValue (maxPower));

  Ssid ssid = Ssid ("AP");
  wifiMac.SetType ("ns3::StaWifiMac",
                   "Ssid", SsidValue (ssid));
  wifiStaDevices.Add (wifi.Install (wifiPhy, wifiMac, wifiStaNodes.Get (0)));

  //Configure the AP node
  wifi.SetRemoteStationManager (manager, "DefaultTxPowerLevel", UintegerValue (powerLevels - 1), "RtsCtsThreshold", UintegerValue (rtsThreshold));
  wifiPhy.Set ("TxPowerStart", DoubleValue (minPower));
  wifiPhy.Set ("TxPowerEnd", DoubleValue (maxPower));
  wifiPhy.Set ("TxPowerLevels", UintegerValue (powerLevels));

  ssid = Ssid ("AP");
  wifiMac.SetType ("ns3::ApWifiMac",
                   "Ssid", SsidValue (ssid));
  wifiApDevices.Add (wifi.Install (wifiPhy, wifiMac, wifiApNodes.Get (0)));
  wifiApDevices.Add (wifi.Install (wifiPhy, wifiMac, wifiApNodes.Get (1)));
  wifiApDevices.Add (wifi.Install (wifiPhy, wifiMac, wifiApNodes.Get (2)));

  wifiDevices.Add (wifiStaDevices);
  wifiDevices.Add (wifiApDevices);

  //Configure the mobility.
  Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator> ();
  //Initial position of AP and STA
  positionAlloc->Add (Vector (ap1_x, ap1_y, 0.0));
  NS_LOG_UNCOND ("Setting initial AP1 position to " << Vector (ap1_x, ap1_y, 0.0));
  positionAlloc->Add (Vector (ap2_x, ap2_y, 0.0));
  NS_LOG_UNCOND ("Setting initial AP2 position to " << Vector (ap2_x, ap2_y, 0.0));
  positionAlloc->Add (Vector (ap3_x, ap3_y, 0.0));
  NS_LOG_UNCOND ("Setting initial AP3 position to " << Vector (ap3_x, ap3_y, 0.0));
  positionAlloc->Add (Vector (sta1_x, sta1_y, 0.0));
  NS_LOG_UNCOND ("Setting initial STA position to " << Vector (sta1_x, sta1_y, 0.0));
  mobility.SetPositionAllocator (positionAlloc);
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobility.Install (wifiApNodes.Get (0));
  mobility.Install (wifiApNodes.Get (1));
  mobility.Install (wifiApNodes.Get (2));
  mobility.Install (wifiStaNodes.Get (0));

  //Configure the IP stack
  InternetStackHelper stack;
  stack.Install (wifiApNodes);
  stack.Install (wifiStaNodes);
  Ipv4AddressHelper address;
  address.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer i = address.Assign (wifiDevices);
  Ipv4Address sinkAddress = i.GetAddress (0);
  uint16_t port = 9;

  //Configure the CBR generator
  PacketSinkHelper sink ("ns3::UdpSocketFactory", InetSocketAddress (sinkAddress, port));
  ApplicationContainer apps_sink = sink.Install (wifiStaNodes.Get (0));

  OnOffHelper onoff ("ns3::UdpSocketFactory", InetSocketAddress (sinkAddress, port));
  onoff.SetConstantRate (DataRate ("54Mb/s"), packetSize);
  onoff.SetAttribute ("StartTime", TimeValue (Seconds (0.5)));
  onoff.SetAttribute ("StopTime", TimeValue (Seconds (simuTime)));
  ApplicationContainer apps_source1 = onoff.Install (wifiApNodes.Get (0));
  ApplicationContainer apps_source2 = onoff.Install (wifiApNodes.Get (1));
  ApplicationContainer apps_source3 = onoff.Install (wifiApNodes.Get (2));

  apps_sink.Start (Seconds (0.1));
  apps_sink.Stop (Seconds (simuTime));

  //Simulator schedule
  Simulator::Schedule (Seconds(0.1), &scheduling, 0.1, 1);

	NS_LOG_UNCOND ("Run Simulation.");
	NS_LOG_UNCOND ("STA moving (1,-0.5) at one time");
  Simulator::Stop (Seconds (simuTime));
  Simulator::Run();
  Simulator::Destroy();
	NS_LOG_UNCOND ("End Simulation.");

  return 0;
}
