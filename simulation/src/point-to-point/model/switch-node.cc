#include "ns3/ipv4.h"
#include "ns3/packet.h"
#include "ns3/ipv4-header.h"
#include "ns3/pause-header.h"
#include "ns3/flow-id-tag.h"
#include "ns3/boolean.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "switch-node.h"
#include "qbb-net-device.h"
#include "ppp-header.h"
#include "ns3/int-header.h"
#include "qbb-header.h"
#include "cn-header.h"
#include "bicc-header.h"
#include <cmath>
#include <ns3/event-id.h>

/**	Flow Table Logging **/
#include <fstream> 
#include <iostream>
#include "ns3/log.h" 
#include "ns3/simulator.h"
NS_LOG_COMPONENT_DEFINE("SwitchNode");
/**	Flow Table Logging **/

/** New L4 Header **/
#include "ns3/custom-header.h"
/** New L4 Header **/

namespace ns3 {

TypeId SwitchNode::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::SwitchNode")
    .SetParent<Node> ()
    .AddConstructor<SwitchNode> ()
	.AddAttribute("EcnEnabled",
			"Enable ECN marking.",
			BooleanValue(false),
			MakeBooleanAccessor(&SwitchNode::m_ecnEnabled),
			MakeBooleanChecker())
	.AddAttribute("CcMode",
			"CC mode.",
			UintegerValue(0),
			MakeUintegerAccessor(&SwitchNode::m_ccMode),
			MakeUintegerChecker<uint32_t>())
	.AddAttribute("AckHighPrio",
			"Set high priority for ACK/NACK or not",
			UintegerValue(0),
			MakeUintegerAccessor(&SwitchNode::m_ackHighPrio),
			MakeUintegerChecker<uint32_t>())
	.AddAttribute("MaxRtt",
			"Max Rtt of the network",
			UintegerValue(9000),
			MakeUintegerAccessor(&SwitchNode::m_maxRtt),
			MakeUintegerChecker<uint32_t>())
	/** New Attribute **/
	.AddAttribute("RPTimer",
			"The rate increase timer at RP in microseconds",
			DoubleValue(200),
			MakeDoubleAccessor(&SwitchNode::t1),
			MakeDoubleChecker<double>())
	.AddAttribute("iBDP",
			"The BDP in the intra-datacenter for BiCC",
			DoubleValue(4000000000),
			MakeDoubleAccessor(&SwitchNode::iBDP),
			MakeDoubleChecker<uint64_t>())
	.AddAttribute("TimeReset",
			"The  reset timer at RP in microseconds",
			DoubleValue(20.0),
			MakeDoubleAccessor(&SwitchNode::TimeReset),
			MakeDoubleChecker<double>())
	.AddAttribute("maxBW",
			"The capacity of bandwitch on DCI switch for BiCC (Mbps)",
			DoubleValue(200000),
			MakeDoubleAccessor(&SwitchNode::maxBW),
			MakeDoubleChecker<uint64_t>())
	.AddAttribute("EwmaGain",
			"Control gain parameter which determines the level of rate decrease",
			DoubleValue(1.0 / 16),
			MakeDoubleAccessor(&SwitchNode::m_g),
			MakeDoubleChecker<double>())
	.AddAttribute ("RateOnFirstCnp",
			"the fraction of rate on first CNP",
			DoubleValue(1.0),
			MakeDoubleAccessor(&SwitchNode::m_rateOnFirstCNP),
			MakeDoubleChecker<double> ())
	.AddAttribute("ClampTargetRate",
			"Clamp target rate.",
			BooleanValue(false),
			MakeBooleanAccessor(&SwitchNode::m_EcnClampTgtRate),
			MakeBooleanChecker())
	.AddAttribute("RPTimerMlx",
			"The rate increase timer at RP in microseconds",
			DoubleValue(1500.0),
			MakeDoubleAccessor(&SwitchNode::m_rpgTimeReset),
			MakeDoubleChecker<double>())
	.AddAttribute("RateDecreaseInterval",
			"The interval of rate decrease check",
			DoubleValue(4.0),
			MakeDoubleAccessor(&SwitchNode::m_rateDecreaseInterval),
			MakeDoubleChecker<double>())
	.AddAttribute("FastRecoveryTimes",
			"The rate increase timer at RP",
			UintegerValue(5),
			MakeUintegerAccessor(&SwitchNode::m_rpgThreshold),
			MakeUintegerChecker<uint32_t>())
	.AddAttribute("AlphaResumInterval",
			"The interval of resuming alpha",
			DoubleValue(55.0),
			MakeDoubleAccessor(&SwitchNode::m_alpha_resume_interval),
			MakeDoubleChecker<double>())
	.AddAttribute("RateAI",
			"Rate increment unit in AI period",
			DataRateValue(DataRate("5Mb/s")),
			MakeDataRateAccessor(&SwitchNode::m_rai),
			MakeDataRateChecker())
	.AddAttribute("RateHAI",
			"Rate increment unit in hyperactive AI period",
			DataRateValue(DataRate("50Mb/s")),
			MakeDataRateAccessor(&SwitchNode::m_rhai),
			MakeDataRateChecker())
	.AddAttribute("MinRate",
			"Minimum rate of a throttled flow",
			DataRateValue(DataRate("100Mb/s")),
			MakeDataRateAccessor(&SwitchNode::m_minRate),
			MakeDataRateChecker())
	.AddAttribute("lineRate",
			"Maximum rate of link",
			DataRateValue(DataRate("100Gb/s")),
			MakeDataRateAccessor(&SwitchNode::lineRate),
			MakeDataRateChecker())
	.AddAttribute("slidingWin",
			"The size of sliding windows for Smooth Start of DCQCN",
			UintegerValue(0),
			MakeUintegerAccessor(&SwitchNode::initWin),
			MakeUintegerChecker<uint32_t>())
	.AddAttribute("initWin",
			"The initial number of tokens",
			UintegerValue(1000),
			MakeUintegerAccessor(&SwitchNode::initWin),
			MakeUintegerChecker<uint32_t>())
	/** DCI_ALG **/
	.AddAttribute("DCI_ALG_FLAG",
		"DCI algorithm flag",
		BooleanValue(false),
		MakeBooleanAccessor(&SwitchNode::m_dciAlgEnabled),
		MakeBooleanChecker())	
	/** DCI_ALG **/		
  /** New Attribute **/
  ;
  
  return tid;
}

SwitchNode::SwitchNode(){
	m_ecmpSeed = m_id;
	m_node_type = 1;
	m_mmu = CreateObject<SwitchMmu>();
	for (uint32_t i = 0; i < pCnt; i++)
		for (uint32_t j = 0; j < pCnt; j++)
			for (uint32_t k = 0; k < qCnt; k++)
				m_bytes[i][j][k] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_txBytes[i] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_lastPktSize[i] = m_lastPktTs[i] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_u[i] = 0;

	/** BICC **/
	// isSendDCI = false;
	// isRecvDCI = false;
	/** BICC **/
	
	forward_table.clear();
	recv_table.clear();

	// through_table.clear();
	through_table2.clear();

	/** BICC **/
	rpTimer = Simulator::Schedule(MicroSeconds(TimeReset), &SwitchNode::CalcEvent, this);	
	counter = 0;	
	/** BICC **/


	/** DCQCN **/
	dcqMap.clear();
	tokenBuckets.clear();
	/** DCQCN **/

	/** Flow Table Logging **/
	/** 初始化 **/
	m_flowTableLoggingEnabled = false;
	m_targetFlowTableSwitchId = 0;
	m_flowTableLogInterval = 1000000.0; // 默认1秒
	m_flowTableLogFilename = "switch_flow_table.log";
	/** Flow Table Logging **/

	/** Control Message **/
	DCI_CM_test = true;
	/** Control Message **/
}

void SwitchNode::CalcEvent()
{
	if (37 == m_mmu->node_id)
	{
		// std::cout << counter << std::endl;
		/*
				if(counter > 0){
					uint64_t val = std::min((uint64_t)(counter * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024), maxBW);
		//			uint64_t val = (uint64_t)(counter * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024);
		//			uint64_t val = counter;
		//			std::cout << clock() << " " << val << std::endl;
					printf("%u %u\n", clock(), val);
				}
				counter = 0;
		*/

		if (through_table2.size() > 0)
		{
			//                      uint64_t val = std::min((uint64_t)(counter * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024), maxBW);
			//                      uint64_t val = (uint64_t)(counter * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024);
			//                      uint64_t val = counter;
			//                      std::cout << clock() << " " << val << std::endl;
			//                      printf("%u %u\n", clock(), val);

			// through_table[ch.sip] += p->GetSize();
			uint64_t totalCnt = 0;
			// std::cout << "====================" << std::endl;
			for (auto kv : through_table2)
			{	
				uint64_t cnt = kv.second;
				uint64_t val = std::min((uint64_t)(cnt * 1.0 * 8 / TimeReset*1000/1024), maxBW);
				totalCnt += val;
				/** Rate Calc **/
				FlowKey key = kv.first;
				m_flowTable[key].currentRate = (uint64_t)(cnt * 1.0 * 8 / TimeReset);
				// std::cout << key.sip << " " << key.dip << " " << key.sport << " " << key.dport << std::endl;
				// std::cout << m_flowTable[key].currentRate << " "<< cnt << std::endl;
				/** Rate Calc **/

				//                                printf("%llu %llu %llu\n", clock(), key, val);
				//  				std::cout << through_table.size() << std::endl;
			}
			std::cout<< "current long haul rate:" << totalCnt<<std::endl;
			// std::cout << "====================" << std::endl;
			/** testDCQCN **
			for(auto kv : dcqMap){
					Mlx mlx = kv.second;
//                              if(mlx.m_targetRate > DataRate(0))
					std::cout << mlx.m_rate.GetBitRate()/1000000000 << " " << mlx.m_targetRate.GetBitRate()/1000000000 << std::endl;
			}

			** testDCQCN **/

			/** testRateLimiting **/
			for (auto kv : dcqMap)
			{
				//                                Mlx mlx = kv.second;
				//                              if(mlx.m_targetRate > DataRate(0))
				// std::cout << mlx.m_rate.GetBitRate()/1000000000 << " " << mlx.m_targetRate.GetBitRate()/1000000000 << std::endl;

				tokenBuckets[kv.first] = std::max((uint64_t)(kv.second.m_targetRate.GetBitRate() * TimeReset / 1000000), m_minRate.GetBitRate());
				//                              std::cout << kv.first << "   " << tokenBuckets[kv.first] << std::endl;
			}
			/** testRateLimiting **/

			// printf("%llu %llu\n", clock(), std::min(totalCnt,maxBW));
			//			std::cout << through_table.size() << std::endl;
			//                        through_table.clear();
		}

		//             counter = 0;
		through_table2.clear();

		Simulator::Cancel(rpTimer);
		rpTimer = Simulator::Schedule(MicroSeconds(TimeReset), &SwitchNode::CalcEvent, this);
	}
}

/** Flow Table **/
/** 用于从自定义包头提取对应 key **/
SwitchNode::FlowKey 
SwitchNode::ExtractFlowKey(CustomHeader &ch) {
    FlowKey key;
    key.sip = ch.sip;
    key.dip = ch.dip;
    key.protocol = ch.l3Prot;
    
    // 根据协议类型设置端口
    if (ch.l3Prot == 0x6) { // TCP
        key.sport = ch.tcp.sport;
        key.dport = ch.tcp.dport;
    } else if (ch.l3Prot == 0x11) { // UDP
        key.sport = ch.udp.sport;
        key.dport = ch.udp.dport;
    } else if (ch.l3Prot == 0xFC || ch.l3Prot == 0xFD) { // ACK/NACK
        key.sport = ch.ack.sport;
        key.dport = ch.ack.dport;
    } else {
        key.sport = 0;
        key.dport = 0;
    }
    
    return key;
}

/** 流表信息打印 **/
void SwitchNode::PrintFlowTable(bool detailed) {
    std::cout << "\n==== Flow Table for Switch " << GetId() << " at " 
              << Simulator::Now().GetSeconds() << "s ====\n";
    std::cout << "Total Active Flows: " << m_flowTable.size() << std::endl;
    
    if (detailed) {
        for (const auto& entry : m_flowTable) {
            const FlowKey& key = entry.first;
            const FlowQueueStats& stats = entry.second;
            
            // 打印流标识信息
            std::cout << "Flow: " 
                      << Ipv4Address(key.sip) << ":" << key.sport << " -> " 
                      << Ipv4Address(key.dip) << ":" << key.dport 
                      << " [Protocol: " << (int)key.protocol << "]\n";
                      
            // 打印流状态信息
            std::cout << "  Queue: " << stats.queueBytes << " bytes, " 
                      << stats.queuePackets << " packets\n";
            std::cout << "  Total: " << stats.totalBytes << " bytes, " 
                      << stats.totalPackets << " packets\n";
            std::cout << "  Last Updated: " << stats.lastUpdateTime.GetSeconds() << "s\n";
            std::cout << "---------------------\n";

			std::cout << "  maxQueuePackets: " << stats.maxQueuePackets << " bytes, " 
                      << stats.maxQueueBytes << " packets\n";
        }
    }
    
    // 计算统计摘要
    uint64_t totalBytes = 0;
    uint64_t totalPackets = 0;
    uint64_t activeQueueBytes = 0;
    uint64_t activeQueuePackets = 0;
    
    for (const auto& entry : m_flowTable) {
        const FlowQueueStats& stats = entry.second;
        totalBytes += stats.totalBytes;
        totalPackets += stats.totalPackets;
        activeQueueBytes += stats.queueBytes;
        activeQueuePackets += stats.queuePackets;
    }
    
    // 打印摘要信息
    std::cout << "Summary Statistics:\n";
    std::cout << "  Active Queue: " << activeQueueBytes << " bytes, " 
              << activeQueuePackets << " packets\n";
    std::cout << "  Cumulative: " << totalBytes << " bytes, " 
              << totalPackets << " packets\n";
    std::cout << "  Avg Bytes Per Flow: " << (m_flowTable.empty() ? 0 : totalBytes / m_flowTable.size()) << "\n";
    std::cout << "===========================\n";
}

/** Flow Table **/


int SwitchNode::GetOutDev(Ptr<const Packet> p, CustomHeader &ch){
//	std::cout << "*****" << std::endl;

	// look up entries
	auto entry = m_rtTable.find(ch.dip);

	// no matching entry
	if (entry == m_rtTable.end())
		return -1;

	// entry found
	auto &nexthops = entry->second;

	// pick one next hop based on hash
	union {
		uint8_t u8[4+4+2+2];
		uint32_t u32[3];
	} buf;
	buf.u32[0] = ch.sip;
	buf.u32[1] = ch.dip;
	if (ch.l3Prot == 0x6)
		buf.u32[2] = ch.tcp.sport | ((uint32_t)ch.tcp.dport << 16);
	else if (ch.l3Prot == 0x11)
		buf.u32[2] = ch.udp.sport | ((uint32_t)ch.udp.dport << 16);
	else if (ch.l3Prot == 0xFC || ch.l3Prot == 0xFD)
		buf.u32[2] = ch.ack.sport | ((uint32_t)ch.ack.dport << 16);

	uint32_t idx = EcmpHash(buf.u8, 12, m_ecmpSeed) % nexthops.size();
	return nexthops[idx];
}

void SwitchNode::CheckAndSendPfc(uint32_t inDev, uint32_t qIndex){
//	if(m_mmu->buffer_size > 15000000) return; // test

	Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
	if (m_mmu->CheckShouldPause(inDev, qIndex)){
		device->SendPfc(qIndex, 0);
		m_mmu->SetPause(inDev, qIndex);
	}

}
void SwitchNode::CheckAndSendResume(uint32_t inDev, uint32_t qIndex){
//	if(m_mmu->buffer_size > 15000000) return; // test
	
	Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
	if (m_mmu->CheckShouldResume(inDev, qIndex)){
		device->SendPfc(qIndex, 1);
		m_mmu->SetResume(inDev, qIndex);
	}
}


/** Control Message **/
void SwitchNode::SendControlMessage(Ptr<Packet> p, uint32_t ifIndex){
	 // 解析原始数据包头部
	 CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
	 ch.getInt = 1; // 解析INT头部
	 p->PeekHeader(ch);

	// 创建控制消息包
    Ptr<Packet> cmPacket = Create<Packet>(0);

	// 创建新的自定义头部
    CustomHeader cmHeader;
    cmHeader.headerType = CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header;

	// 设置控制消息头部内容
    cmHeader.DCI_CM_header.test = 123;  // 设置测试字段

	/** IPV4_H **/	
	// 交换源目的IP地址
    cmHeader.sip = ch.dip;
    cmHeader.dip = ch.sip;

	// 设置为DCI控制消息类型
    cmHeader.l3Prot = DCI_CM;  // 0xFA

	// 设置IPv4头部其他字段
    cmHeader.m_ttl = 64;
    cmHeader.ipv4Flags = 0;
    cmHeader.m_payloadSize = 0;  // 没有payload
	/** IPV4_H **/

	// PPP协议号
	cmHeader.pppProto = 0x0021;  

	// 准备发送
	cmPacket->AddHeader(cmHeader);
    cmPacket->AddPacketTag(FlowIdTag(ifIndex));

	// 发送控制消息
	int idx = GetOutDev(cmPacket, cmHeader);
	std::cout << idx << "test " << std::endl;
	m_devices[idx]->SwitchSend(0, cmPacket, cmHeader, true);

}
/** Control Message **/

/** SRC-DCI CNP GENERATION **/
// 生成DCI CNP报文
void SwitchNode::SrcDCICNPGen(Ptr<Packet> p, uint32_t ifIndex) {
	// 解析原始数据包头部
	CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
	ch.getInt = 1; // 解析INT头部
	p->PeekHeader(ch);

   // 创建控制消息包
   Ptr<Packet> cmPacket = Create<Packet>(0);

   // 创建新的自定义头部
   CustomHeader cmHeader;
   cmHeader.headerType = CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header;

   // 设置控制消息头部内容
   cmHeader.SRC_DCI_CNP_header.sport = ch.udp.dport; 
   cmHeader.SRC_DCI_CNP_header.dport = ch.udp.sport;
   cmHeader.SRC_DCI_CNP_header.pg = ch.udp.pg;

   /** IPV4_H **/	
   // 交换源目的IP地址
   cmHeader.sip = ch.dip;
   cmHeader.dip = ch.sip;

   // 设置为DCI控制消息类型
   cmHeader.l3Prot = SRC_DCI_CNP;  

   // 设置IPv4头部其他字段
   cmHeader.m_ttl = 64;
   cmHeader.ipv4Flags = 0;
   cmHeader.m_payloadSize = 0;  // 没有payload
   /** IPV4_H **/

   // PPP协议号
   cmHeader.pppProto = 0x0021;  

   // 准备发送
   cmPacket->AddHeader(cmHeader);
   cmPacket->AddPacketTag(FlowIdTag(ifIndex));

   // 发送控制消息
   int idx = GetOutDev(cmPacket, cmHeader);
//    std::cout << "SRC-DCI CNP GENERATION " << m_mmu->node_id << std::endl;
   m_devices[idx]->SwitchSend(0, cmPacket, cmHeader, true);

}
/** SRC-DCI CNP GENERATION **/


void SwitchNode::SendToDev(Ptr<Packet>p, CustomHeader &ch){
	/**	Flow Table**/
	/** 包到达处对对应流表项更新 **/
	if (ch.l3Prot == 0x11 && (m_mmu->node_id == DCI_SWITCH_0 || m_mmu->node_id == DCI_SWITCH_1)){
		 // 在 DCI 交换机处针对UDP流量进行信息统计
		FlowKey key= ExtractFlowKey(ch);
		// 查找或创建流记录
		auto it = m_flowTable.find(key);
		if (it == m_flowTable.end()) {
			// 新流
			FlowQueueStats stats;
			stats.queueBytes = p->GetSize();
			stats.queuePackets = 1;
			stats.totalBytes = p->GetSize();
			stats.totalPackets = 1;
			stats.lastUpdateTime = Simulator::Now();
			m_flowTable[key] = stats;
			stats.maxQueueBytes = p->GetSize();  // 初始值即为当前值
        	stats.maxQueuePackets = 1;           // 初始值即为当前值
			stats.lastCnpTime = Seconds(-1.0);   // 初始化为负
		} else {
			// 更新现有流
			it->second.queueBytes += p->GetSize();
			it->second.queuePackets++;
			it->second.totalBytes += p->GetSize();
			it->second.totalPackets++;
			it->second.lastUpdateTime = Simulator::Now();

			// 更新最大队列大小
			if (it->second.queueBytes > it->second.maxQueueBytes) {
				it->second.maxQueueBytes = it->second.queueBytes;
			}
			if (it->second.queuePackets > it->second.maxQueuePackets) {
				it->second.maxQueuePackets = it->second.queuePackets;
			}
		}
	}
	/**	Flow Table**/

	/** Control Message **/
	// 经测试可以在 0 号 DCI 交换机上收到控制报文
	// if (ch.l3Prot ==0xFA){
	// 	std::cout << "Switch" <<GetId() << std::endl;
	// 	std::cout << "DCI CM" << (int)ch.DCI_CM_header.test<<std::endl;
	// 	return; // Drop
	// }
	/** Control Message **/

	int idx = GetOutDev(p, ch);

	/** Through Table **/
	// 按五元组统计不同流的流量
	// std::cout << "*****" << std::endl;
	// counter += p->GetSize();
	
	// std::cout << ch.sip << std::endl;
	// through_table[ch.sip] += p->GetSize();
	// std::cout << ch.sip  << std::endl;
	// if (ch.l3Prot == 0x11 || ch.l3Prot == 0x6)
	// {
	// 	// std::cout << ch.l3Prot << " "  << ch.sip << std::endl;
	// 	// through_table[ch.sip] += p->GetSize();
	// 	// uint64_t flowId = ((uint64_t)ch.sip << 32) | ((uint64_t)ch.udp.pg << 16) | (uint64_t)ch.udp.sport;
	// 	FlowKey key = ExtractFlowKey(ch);

	// 	// std::cout << key << std::endl;
	// 	// through_table[flowId] += p->GetSize();
	// 	through_table2[key] += p->GetSize();
	// }
	// else if(ch.l3Prot == 0xFC || ch.l3Prot == 0xFD){ // || ch.l3Prot == 0xFF || ch.l3Prot == 0xFE
    // 	// std::cout << ch.l3Prot << " "  << ch.dip << std::endl;
    //     // through_table[ch.dip] += p->GetSize();
    //     // uint64_t flowId = ((uint64_t)ch.dip << 32) | ((uint64_t)ch.ack.pg << 16) | (uint64_t)ch.ack.dport;
    //     FlowKey key= ExtractFlowKey(ch);
		
	// 	// std::cout << key << std::endl;
    //     // through_table[flowId] += p->GetSize();
	// 	through_table2[key] += p->GetSize();
    // }
	/** Through Table **/


	if (idx >= 0){
		NS_ASSERT_MSG(m_devices[idx]->IsLinkUp(), "The routing table look up should return link that is up");

		// determine the qIndex
		uint32_t qIndex;
		if (ch.l3Prot == 0xFF || ch.l3Prot == 0xFE || (m_ackHighPrio && (ch.l3Prot == 0xFD || ch.l3Prot == 0xFC))){  //QCN or PFC or NACK, go highest priority
			qIndex = 0;
		}else{
			qIndex = (ch.l3Prot == 0x06 ? 1 : ch.udp.pg); // if TCP, put to queue 1
		}

		// admission control
		FlowIdTag t;
		p->PeekPacketTag(t);
		uint32_t inDev = t.GetFlowId();
		if (qIndex != 0){ //not highest priority
			if (m_mmu->CheckIngressAdmission(inDev, qIndex, p->GetSize()) && m_mmu->CheckEgressAdmission(idx, qIndex, p->GetSize())){			// Admission control
				m_mmu->UpdateIngressAdmission(inDev, qIndex, p->GetSize());
				m_mmu->UpdateEgressAdmission(idx, qIndex, p->GetSize());
			}else{
				return; // Drop
			}
			/** TEST **/
			/** test **
			if(m_mmu->buffer_size > 15000000 && (false == m_mmu->paused[inDev][qIndex]) 
					&& (m_mmu->GetSharedUsed(inDev, qIndex) > 0)){
				Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
				if(false == m_mmu->paused[inDev][qIndex]){
					m_mmu->SetPause(inDev, qIndex);
					Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                	device->SendPfc(qIndex, 0);
				}

			//	std::cout << "Pause:  " << inDev << " " << qIndex << " " 
			//		<< m_mmu->GetSharedUsed(inDev, qIndex) << std::endl;
			}

	                if(m_mmu->buffer_size > 15000000 && (true == m_mmu->paused[inDev][qIndex])
                                && (m_mmu->GetSharedUsed(inDev, qIndex) < m_mmu->pfc_remote)){
        	                Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                	        if(true == m_mmu->paused[inDev][qIndex]){
                        	        m_mmu->SetResume(inDev, qIndex);
                                	Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                	device->SendPfc(qIndex, 1);
                        	}

                       		// std::cout << "Resume:  " << inDev << " " << qIndex << " "
                        	//      << m_mmu->GetSharedUsed(inDev, qIndex) << std::endl;
                	}

			** test **/
			/** test **
			if(m_mmu->buffer_size > 15000000){
				if(false == m_mmu->paused[inDev][qIndex]){
					if(m_mmu->GetSharedUsed(inDev, qIndex) > 0){
						Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
		                                if(false == m_mmu->paused[inDev][qIndex]){
                		                        m_mmu->SetPause(inDev, qIndex);
                                		        Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                        		device->SendPfc(qIndex, 0);
                                		}
        	        	 //       	      std::cout << "Pause:  " << inDev << " " << qIndex << " "
	        	         //	                      << m_mmu->GetSharedUsed(inDev, qIndex) << std::endl;
					}
				}
				else{
					if(m_mmu->GetSharedUsed(inDev, qIndex) < m_mmu->pfc_remote * 20){ // m_mmu->pfc_remote
						Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                		if(true == m_mmu->paused[inDev][qIndex]){
                                        		m_mmu->SetResume(inDev, qIndex);
                                        		Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                        		device->SendPfc(qIndex, 1);
                                		}

                               // 		 std::cout << "Resume:  " << inDev << " " << qIndex << " "
                               // 		      << m_mmu->GetSharedUsed(inDev, qIndex) << std::endl;
                                        }
				}
			}

			** test **/
			/*
			if(m_mmu->buffer_size > 15000000 && m_mmu->GetSharedUsed(inDev, qIndex)>0 
					&& false == m_mmu->paused[inDev][qIndex]){

				// m_mmu->SetPause(inDev, qIndex);
                                Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                
				if(false == m_mmu->paused[inDev][qIndex]){
					device->SendPfc(qIndex, 0);
				}

				m_mmu->SetPause(inDev, qIndex);
			}
			*/
		
		/*
			if(m_mmu->buffer_size > 15000000 && m_mmu->GetSharedUsed(inDev, qIndex)>0)
				std::cout <<  inDev << " " << qIndex << " " << m_mmu->GetSharedUsed(inDev, qIndex) << std::endl;
		*/		
		/** TEST **/	
			CheckAndSendPfc(inDev, qIndex);
		}

		/* BICC */
		uint8_t ecnbits = ch.GetIpv4EcnBits();
        // bool egressCongested = m_ecnEnabled && m_mmu->ShouldSendCN(idx, qIndex);
        // if((m_dciAlgEnabled ||9 == m_ccMode) && (DCI_SWITCH_0== m_mmu->node_id || DCI_SWITCH_1== m_mmu->node_id) && ch.l3Prot == 0x11 && (egressCongested || ecnbits)){
		/** NEW ECN LOGIC **/
		// 解释：这里只会在 DCI 交换机上为源 DC 内已有 ECN 标记的报文进行 CNP 代理，DCI 交换机本身的 ECN 标记 代理在 SwitchNotifyDequeue 中实现
		if((m_dciAlgEnabled ||9 == m_ccMode) && (DCI_SWITCH_0== m_mmu->node_id || DCI_SWITCH_1== m_mmu->node_id) && ch.l3Prot == 0x11 && ecnbits){
		/** NEW ECN LOGIC **/
			/** CNP SELECT **/
			if (m_ccMode == 9)
			{
				sendCNPByDCI(p, idx);
			}
			else if (m_ccMode == 1)
			{
				SrcDCICNPGen(p, idx);
			}
			/** CNP SELECT **/
			std::cout << "ingress SRC-DCI CNP ECN Clean" << m_mmu->node_id <<std::endl;

			PppHeader ppp;
            Ipv4Header h;
            p->RemoveHeader(ppp);
            p->RemoveHeader(h);
            h.SetEcn((Ipv4Header::EcnType)0x00);
            p->AddHeader(h);
            p->AddHeader(ppp);
		}
		/* BICC */

		/** Control Message **/
		// 先测试报文构造可行性
		// if (DCI_CM_test && m_mmu->node_id == 37)
		// {
		// 	// FlowKey key= ExtractFlowKey(ch);
		// 	// SendControlMessage(p,idx);
		// 	DCI_CM_test = false;
		// 	// std::cout << "DCI CM" << std::endl;
		// 	// std::cout << "DCI_FLAG "<<m_dciAlgEnabled << std::endl;
		// 	std::cout <<" DCI_BUFFER_SIZE "<< m_mmu->buffer_size << std::endl;
		// 	std::cout <<" DCI RESERVE "<< m_mmu->reserve << std::endl;
		// 	std::cout <<" DCI resume_offset "<< m_mmu->resume_offset << std::endl;
		// 	// std::cout << "mccmode"<< m_ccMode << std::endl;
		// }
		/** Control Message **/

		/** SRC-DCI CNP TEST **/
		// 测试CNP报文有效性
		// if (ch.l3Prot == 0x11 && m_mmu->node_id == DCI_SWITCH_1)
		// {
		// 	SrcDCICNPGen(p,idx);
		// }
		/** SRC-DCI CNP TEST **/

		/* BICC */
		/*
		if(1 == m_mmu->node_id && ch.l3Prot == 0x11 && ecnbits){
			PppHeader ppp;
			Ipv4Header h;
			p->RemoveHeader(ppp);
			p->RemoveHeader(h);
			h.SetEcn((Ipv4Header::EcnType)0x00);
			p->AddHeader(h);
			p->AddHeader(ppp);

			// std::cout << "111111111111" << std::endl;
		}
		*/
		/* BICC */

		m_bytes[inDev][idx][qIndex] += p->GetSize();

		/** TEST **/
		// counter += p->GetSize();		
		// m_devices[idx]->SwitchSend(qIndex, p, ch);
		/** TEST **/

		/** BiCC third loop **/
		Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(m_devices[idx]);
		uint32_t fip = ch.dip;
		bool isSend = true;

		if(9 == m_ccMode && 37 == m_mmu->node_id && ch.l3Prot == 0x11){
			/** Remove ECN Signnal **/
			PppHeader ppp;
            Ipv4Header h;
            p->RemoveHeader(ppp);
            p->RemoveHeader(h);
            h.SetEcn((Ipv4Header::EcnType)0x00);
            p->AddHeader(h);
            p->AddHeader(ppp);
			/** Remove ECN Signnal **/

            /** token Rate Limiting **/
            // uint64_t flowId = ((uint64_t)ch.sip << 32) | ((uint64_t)ch.udp.pg << 16) | (uint64_t)ch.udp.sport;
            FlowKey flowId = ExtractFlowKey(ch); 
			
			if(tokenBuckets.find(flowId) != tokenBuckets.end()){
                if(tokenBuckets[flowId] >= p->GetSize()){
                    tokenBuckets[flowId] -= p->GetSize();
                }else{
					 std::cout << tokenBuckets[flowId] << " ******1******"  << std::endl;
 					/*
                    for(auto kv : dcqMap){
                        std::cout << kv.first << "   " << kv.second.m_targetRate << "  " << kv.second.m_rate  << std::endl;
                    }
					*/
                    isSend = false;
                }
            }
        	/** token Rate Limiting **/

            if(forward_table.find(fip)!=forward_table.end()
                && recv_table.find(fip)!=recv_table.end()
                && forward_table[fip]>iBDP){
            	// std::cout << recv_table[fip] << " ******1******" << std::endl;

                isSend = false;
            }else{
                forward_table[fip] += p->GetSize();
            }


			/*
			if(forward_table.find(fip)!=forward_table.end() 
					&& recv_table.find(fip)!=recv_table.end()
				      	&& forward_table[fip]>iBDP){
				//  std::cout << recv_table[fip] << " ******1******" << std::endl;
				return;
			}

			forward_table[fip] += p->GetSize();
			// std::cout << bi_table[fip] << " ******1******" << std::endl;
			*/
		}

		if (9 == m_ccMode && 37 == m_mmu->node_id && ch.l3Prot == 0xFC)
		{
			// std::cout << ch.udp.ih.recvBiCCBytes << " ******1******" << std::endl;
			// uint32_t aSip = ch.sip;
			recv_table[ch.sip] = ch.ack.ih.recvBiCCBytes;
			// if(bi_table[aSip] < ch.udp.ih.recvBiCCBytes)
			// std::cout << forward_table[aSip]  << " ******2******" << recv_table[aSip]  << std::endl;
			// if(bi_table.find(fip) != bi_table.end()){
			// 	bi_table[fip] -= ch.udp.ih.recvBiCCBytes;
			// 	std::cout << bi_table[fip] << std::endl;

			// 	if(bi_table[fip] <= 0){
			// 		bi_table.erase(fip);
			// 	}
			// }
			// ch.udp.ih.recvBiCCBytes = 1;
			// std::cout << ch.udp.ih.isLongLoop << " ******2******" << std::endl;

			/** testDCQCN **/
			// uint64_t flowId = ((uint64_t)ch.dip << 32) | ((uint64_t)ch.ack.pg << 16) | (uint64_t)ch.ack.dport;
			FlowKey flowId = ExtractFlowKey(ch);
			if (dcqMap.find(flowId) == dcqMap.end())
			{
				Mlx mlx;
				mlx.m_first_cnp = true;
				mlx.m_rate = lineRate;
				mlx.m_targetRate = lineRate;

				dcqMap[flowId] = mlx;
			}

			// uint8_t cnp = (ch.ack.flags >> qbbHeader::FLAG_CNP) & 1;

			if (ch.ack.flags == 1)
			{
				// std::cout << ch.ack.flags << std::endl;
				// uint64_t flowId = ((uint64_t)ch.dip << 32) | ((uint64_t)ch.ack.pg << 16) | (uint64_t)ch.ack.dport;
				if (0 == slidingWin)
				{
					cnp_received_mlx(flowId);
					slidingWin = initWin;
				}
				else
				{
					slidingWin -= 1;
				}
			}
			/** testDCQCN **/
		}
		/** BiCC third loop **/	

		/** New Send **/
		// 通过 isSend 标记确认是否需要发送
		/*
			if(isSend) {
					m_devices[idx]->SwitchSend(qIndex, p, ch, true);
				}
		)
		*/
		// m_devices[idx]->SwitchSend(qIndex, p, ch, isSend);
		m_devices[idx]->SwitchSend(qIndex, p, ch, true);// 测试用，暂不考虑数据包缓存，直接转发

		/** New Send **/

		/** Original Send **/
		// m_devices[idx]->SwitchSend(qIndex, p, ch);
		/** Original Send **/

		/** TEST **/
		// std::cout << dev->GetDataRate().GetBitRate() << std::endl;
		/** TEST **/
	}else
		return; // Drop
}


uint32_t SwitchNode::EcmpHash(const uint8_t* key, size_t len, uint32_t seed) {
  uint32_t h = seed;
  if (len > 3) {
    const uint32_t* key_x4 = (const uint32_t*) key;
    size_t i = len >> 2;
    do {
      uint32_t k = *key_x4++;
      k *= 0xcc9e2d51;
      k = (k << 15) | (k >> 17);
      k *= 0x1b873593;
      h ^= k;
      h = (h << 13) | (h >> 19);
      h += (h << 2) + 0xe6546b64;
    } while (--i);
    key = (const uint8_t*) key_x4;
  }
  if (len & 3) {
    size_t i = len & 3;
    uint32_t k = 0;
    key = &key[i - 1];
    do {
      k <<= 8;
      k |= *key--;
    } while (--i);
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    h ^= k;
  }
  h ^= len;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

void SwitchNode::SetEcmpSeed(uint32_t seed){
	m_ecmpSeed = seed;
}

void SwitchNode::AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx){
	uint32_t dip = dstAddr.Get();
	m_rtTable[dip].push_back(intf_idx);
}

void SwitchNode::ClearTable(){
	m_rtTable.clear();
}

// This function can only be called in switch mode
bool SwitchNode::SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch){
	SendToDev(packet, ch);

	/** BICC **/
        // std::cout << m_mmu->node_id << std::endl;
        /* if(0 == m_mmu->node_id && ch.l3Prot == 0x11 && ch.udp.ih.rate != 0){
		* Send CNP *
                qbbHeader seqh;
                seqh.SetSeq(0);
                seqh.SetPG(0);
                seqh.SetSport(ch.udp.dport);
                seqh.SetDport(ch.udp.sport);

                IntHeader ih;
                seqh.SetIntHeader(ih);

                seqh.SetStart(6);

                Ptr<Packet> CNP_Pkt = Create<Packet>(std::max(60-14-20-(int)seqh.GetSerializedSize(), 0));
                CNP_Pkt->AddHeader(seqh);

                Ipv4Header head;        // Prepare IPv4 header
                head.SetDestination(Ipv4Address(ch.sip));
                head.SetSource(Ipv4Address(ch.dip));
                head.SetProtocol(0xFB); //cnp
                head.SetTtl(64);
                head.SetPayloadSize(CNP_Pkt->GetSize());
                head.SetIdentification(0);
                CNP_Pkt->AddHeader(head);

		PppHeader ppp;
		ppp.SetProtocol(0x0021);
		CNP_Pkt->AddHeader(ppp);

		CustomHeader chh(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
		CNP_Pkt->PeekHeader(chh);
		CNP_Pkt->AddPacketTag(FlowIdTag(0));

                // send
	        SendToDev(CNP_Pkt, chh);               
	}*/
        /** BICC **/

	return true;
}

void SwitchNode::SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p){
	FlowIdTag t;
	p->PeekPacketTag(t);
	counter += p->GetSize(); 
	
	/** Through Table **/
	// through_table[ch.sip] += p->GetSize();	
	// uint32_t sip = 0x0b000001 + ((m_id / 256) * 0x00010000) + ((m_id % 256) * 0x00000100);
	CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
    ch.getInt = 1; // parse INT header
    p->PeekHeader(ch);
	// through_table[ch.sip] += p->GetSize();
	
	if (ch.l3Prot == 0x11 || ch.l3Prot == 0x6){
		// std::cout << ch.l3Prot << " "  << ch.sip << std::endl;	
		// through_table[ch.sip] += p->GetSize();
		// uint64_t flowId = ((uint64_t)ch.sip << 32) | ((uint64_t)ch.udp.pg << 16) | (uint64_t)ch.udp.sport;
        FlowKey key= ExtractFlowKey(ch);

    	// std::cout << key << std::endl;
    	// through_table[flowId] += p->GetSize();
    	through_table2[key] += p->GetSize();		
		// through_table[flowId] += p->GetSize();
	}
	// else if(ch.l3Prot == 0xFC || ch.l3Prot == 0xFD){ // || ch.l3Prot == 0xFF || ch.l3Prot == 0xFE
	// 	// std::cout << ch.l3Prot << " "  << ch.dip << std::endl;
	// 	// through_table[ch.dip] += p->GetSize();
    // 	// uint64_t flowId = ((uint64_t)ch.dip << 32) | ((uint64_t)ch.ack.pg << 16) | (uint64_t)ch.ack.dport;
    // 	// through_table[flowId] += p->GetSize();
    //     FlowKey key= ExtractFlowKey(ch);

    //     // std::cout << key << std::endl;
    //     // through_table[flowId] += p->GetSize();
    //     through_table2[key] += p->GetSize();
	// }
	/** Through Table **/

	/** Flow Table **/
	/** 包发送后对对应流表项更新 **/
	if (ch.l3Prot == 0x11){
		FlowKey key= ExtractFlowKey(ch);
		// 查找流记录
		auto it = m_flowTable.find(key);
		if (it != m_flowTable.end()) {
			// 更新现有流
			it->second.queueBytes -= p->GetSize();
			it->second.queuePackets--;
			it->second.lastUpdateTime = Simulator::Now();
			// if (it->second.queueBytes == 0) {
			// 	m_flowTable.erase(it); // 删除流记录
			// }
		}
	}
	/** Flow Table **/

	if (qIndex != 0){
		uint32_t inDev = t.GetFlowId();
		m_mmu->RemoveFromIngressAdmission(inDev, qIndex, p->GetSize());
		m_mmu->RemoveFromEgressAdmission(ifIndex, qIndex, p->GetSize());
		m_bytes[inDev][ifIndex][qIndex] -= p->GetSize();
		if (m_ecnEnabled){
			bool egressCongested = m_mmu->ShouldSendCN(ifIndex, qIndex);
			/** TEST **/
			//std::cout << "22222222222" << std::endl;
			/** TEST **/

			if (egressCongested)
			{
				/** Bug Fix **/
				PppHeader ppp;
				Ipv4Header h;
				p->RemoveHeader(ppp);
				p->RemoveHeader(h);
				h.SetEcn((Ipv4Header::EcnType)0x03);
				p->AddHeader(h);
				p->AddHeader(ppp);
				/** TEST **/
				std::cout << "Congestion Marked At Switch" << m_mmu->node_id <<std::endl;
				/** TEST **/
				/** Bug Fix **/
				// DCI Egress 处不清除 ECN 标记 -> ECN 作为一种拥塞标志是为了概率减少交换机队列长度
				// if ((m_mmu->node_id == DCI_SWITCH_0 || m_mmu->node_id == DCI_SWITCH_1) && (m_ccMode == 9||m_dciAlgEnabled) && ch.l3Prot == 0x11)
				// {
				// 	if (m_dciAlgEnabled && m_ccMode == 1)
				// 	{
				// 		SrcDCICNPGen(p, ifIndex);
				// 		// std::cout << "SRC-DCI CNP Egress Clean" << std::endl;
				// 	}
				// 	else if (m_ccMode == 9)
				// 	{
				// 		sendCNPByDCI(p, ifIndex);
				// 	}
				// 	PppHeader ppp;
				// 	Ipv4Header h;
				// 	p->RemoveHeader(ppp);
				// 	p->RemoveHeader(h);
				// 	h.SetEcn((Ipv4Header::EcnType)0x00);
				// 	p->AddHeader(h);
				// 	p->AddHeader(ppp);
				// }

				/** BICC **/
				// std::cout << m_mmu->node_id << std::endl;
				// if(0 == m_mmu->node_id){
					// 
				// }
				/** BICC **/
			}
		}       
		/** TEST **/
		/*
		if(m_mmu->buffer_size > 15000000 && m_mmu->paused[inDev][qIndex]){
			std::cout << inDev << " " << qIndex << " "  << m_mmu->GetSharedUsed(inDev, qIndex) << " " << m_mmu->pfc_remote * 0.01 << std::endl;
		}
		*/
		/*                
		if(m_mmu->buffer_size > 15000000 && (true == m_mmu->paused[inDev][qIndex])
                                && m_mmu->GetSharedUsed(inDev, qIndex) < m_mmu->pfc_remote){
			Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                        if(true == m_mmu->paused[inDev][qIndex]){
				m_mmu->SetResume(inDev, qIndex);
                                Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                device->SendPfc(qIndex, 1);
                        }

                       // std::cout << "Resume:  " << inDev << " " << qIndex << " " 
			//	<< m_mmu->GetSharedUsed(inDev, qIndex) << std::endl;
                }
		*/
                /** test  1047754 **/
		/*
		if(m_mmu->buffer_size > 15000000 && true == m_mmu->paused[inDev][qIndex] 
				&& m_mmu->GetSharedUsed(inDev, qIndex) <= 104775){

			// if(lastBytes>0 && lastBytes>m_mmu->GetSharedUsed(inDev, qIndex)){
				// std::cout << inDev << " " << qIndex << " "  << m_mmu->paused[inDev][qIndex] << " " << m_mmu->GetSharedUsed(inDev, qIndex) << " " << m_mmu->hdrm_bytes[inDev][qIndex] << std::endl;
                                Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                
				if(true == m_mmu->paused[inDev][qIndex]){
					device->SendPfc(qIndex, 1);	
				}

				m_mmu->SetResume(inDev, qIndex);
			// }
			
			// lastBytes = m_mmu->GetSharedUsed(inDev, qIndex);
		}
		*/
		/** TEST **/	
		//CheckAndSendPfc(inDev, qIndex);
		CheckAndSendResume(inDev, qIndex);
	}
	if (1){
		uint8_t* buf = p->GetBuffer();
		if (buf[PppHeader::GetStaticSize() + 9] == 0x11){ // udp packet
			IntHeader *ih = (IntHeader*)&buf[PppHeader::GetStaticSize() + 20 + 8 + 6]; // ppp, ip, udp, SeqTs, INT
			Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(m_devices[ifIndex]);
			if (m_ccMode == 3 || m_ccMode == 1 || m_ccMode == 7 || m_ccMode == 9){ // HPCC DCQCN ** test **
				ih->PushHop(Simulator::Now().GetTimeStep(), m_txBytes[ifIndex], dev->GetQueue()->GetNBytesTotal(), dev->GetDataRate().GetBitRate());
			}else if (m_ccMode == 10){ // HPCC-PINT
				uint64_t t = Simulator::Now().GetTimeStep();
				uint64_t dt = t - m_lastPktTs[ifIndex];
				if (dt > m_maxRtt)
					dt = m_maxRtt;
				uint64_t B = dev->GetDataRate().GetBitRate() / 8; //Bps
				uint64_t qlen = dev->GetQueue()->GetNBytesTotal();
				double newU;

				/**************************
				 * approximate calc
				 *************************/
				int b = 20, m = 16, l = 20; // see log2apprx's paremeters
				int sft = logres_shift(b,l);
				double fct = 1<<sft; // (multiplication factor corresponding to sft)
				double log_T = log2(m_maxRtt)*fct; // log2(T)*fct
				double log_B = log2(B)*fct; // log2(B)*fct
				double log_1e9 = log2(1e9)*fct; // log2(1e9)*fct
				double qterm = 0;
				double byteTerm = 0;
				double uTerm = 0;
				if ((qlen >> 8) > 0){
					int log_dt = log2apprx(dt, b, m, l); // ~log2(dt)*fct
					int log_qlen = log2apprx(qlen >> 8, b, m, l); // ~log2(qlen / 256)*fct
					qterm = pow(2, (
								log_dt + log_qlen + log_1e9 - log_B - 2*log_T
								)/fct
							) * 256;
					// 2^((log2(dt)*fct+log2(qlen/256)*fct+log2(1e9)*fct-log2(B)*fct-2*log2(T)*fct)/fct)*256 ~= dt*qlen*1e9/(B*T^2)
				}
				if (m_lastPktSize[ifIndex] > 0){
					int byte = m_lastPktSize[ifIndex];
					int log_byte = log2apprx(byte, b, m, l);
					byteTerm = pow(2, (
								log_byte + log_1e9 - log_B - log_T
								)/fct
							);
					// 2^((log2(byte)*fct+log2(1e9)*fct-log2(B)*fct-log2(T)*fct)/fct) ~= byte*1e9 / (B*T)
				}
				if (m_maxRtt > dt && m_u[ifIndex] > 0){
					int log_T_dt = log2apprx(m_maxRtt - dt, b, m, l); // ~log2(T-dt)*fct
					int log_u = log2apprx(int(round(m_u[ifIndex] * 8192)), b, m, l); // ~log2(u*512)*fct
					uTerm = pow(2, (
								log_T_dt + log_u - log_T
								)/fct
							) / 8192;
					// 2^((log2(T-dt)*fct+log2(u*512)*fct-log2(T)*fct)/fct)/512 = (T-dt)*u/T
				}
				newU = qterm+byteTerm+uTerm;

				#if 0
				/**************************
				 * accurate calc
				 *************************/
				double weight_ewma = double(dt) / m_maxRtt;
				double u;
				if (m_lastPktSize[ifIndex] == 0)
					u = 0;
				else{
					double txRate = m_lastPktSize[ifIndex] / double(dt); // B/ns
					u = (qlen / m_maxRtt + txRate) * 1e9 / B;
				}
				newU = m_u[ifIndex] * (1 - weight_ewma) + u * weight_ewma;
				printf(" %lf\n", newU);
				#endif

				/************************
				 * update PINT header
				 ***********************/
				uint16_t power = Pint::encode_u(newU);
				if (power > ih->GetPower())
					ih->SetPower(power);

				m_u[ifIndex] = newU;
			}
		}
	}
	m_txBytes[ifIndex] += p->GetSize();
	m_lastPktSize[ifIndex] = p->GetSize();
	m_lastPktTs[ifIndex] = Simulator::Now().GetTimeStep();
}

int SwitchNode::logres_shift(int b, int l){
	static int data[] = {0,0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
	return l - data[b];
}

int SwitchNode::log2apprx(int x, int b, int m, int l){
	int x0 = x;
	int msb = int(log2(x)) + 1;
	if (msb > m){
		x = (x >> (msb - m) << (msb - m));
		#if 0
		x += + (1 << (msb - m - 1));
		#else
		int mask = (1 << (msb-m)) - 1;
		if ((x0 & mask) > (rand() & mask))
			x += 1<<(msb-m);
		#endif
	}
	return int(log2(x) * (1<<logres_shift(b, l)));
}

void SwitchNode::sendCNPByDCI(Ptr<Packet> p, uint32_t ifIndex){
//	std::cout << "SendDCI send"<< std::endl;
	/** test **/
//	if(1) return;	

	CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
	ch.getInt = 1; // parse INT header
	p->PeekHeader(ch);
	Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(m_devices[ifIndex]);
	Ptr<Packet> newp = Create<Packet>(0);
/*
	qbbHeader seqh;
	seqh.SetSeq(0);
	seqh.SetPG(ch.udp.pg);
	seqh.SetSport(ch.udp.dport);
	seqh.SetDport(ch.udp.sport);
	seqh.SetIntHeader(ch.udp.ih);		
	newp->AddHeader(seqh);
*/
	BICCHeader bicc;
	bicc.qIndex = ch.udp.pg;
	bicc.dport = ch.udp.sport;
	newp->AddHeader(bicc);
	
/*	
	CnHeader cnp;
	cnp.SetQindex(ch.udp.sport);
	cnp.SetFlow(ch.udp.sport);
	cnp.SetECNBits(1);
	cnp.SetQfb(1100);
	cnp.SetTotal(1100);
*/

//	bicc.dport = ch.udp.sport;
//	newp->AddHeader(bicc);

	Ipv4Header ipv4h; // Prepare IPv4 header
	ipv4h.SetProtocol(0xFB);
	uint32_t sip = 0x0b000001 + ((m_id / 256) * 0x00010000) + ((m_id % 256) * 0x00000100);
	ipv4h.SetSource(Ipv4Address(ch.dip));
	ipv4h.SetDestination(Ipv4Address(ch.sip));
//	std::cout << ch.udp.sport << std::endl;
	Ipv4Address dip = Ipv4Address(ch.sip);
	ipv4h.SetPayloadSize(0);
	ipv4h.SetTtl(64);
	newp->AddHeader(ipv4h);
	PppHeader ppp;
	ppp.SetProtocol(0x0021);
	newp->AddHeader(ppp);
	CustomHeader chh(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
	newp->PeekHeader(chh);
	newp->AddPacketTag(FlowIdTag(ifIndex));
	
	SendToDev(newp, chh);
}

/** Flow Table Logging**/
/** 用于检测不同流队列长度等信息是否统计正确 -- 已验证**/
/** 初始化 Logging **/
void SwitchNode::StartFlowTableLogging(uint32_t targetSwitchId, double intervalMs, std::string filename) {
    m_flowTableLoggingEnabled = true;
    m_targetFlowTableSwitchId = targetSwitchId;
    m_flowTableLogFilename = filename;
    m_flowTableLogInterval = intervalMs * 1000.0; // 转换为微秒
    
    // 立即执行一次记录，并安排定期执行
    LogFlowTablePeriodically();
}

/** Logging 逻辑 **/
void SwitchNode::LogFlowTablePeriodically() {
    // 如果不是目标交换机或功能被禁用，则直接返回
    if (!m_flowTableLoggingEnabled || m_mmu->node_id != m_targetFlowTableSwitchId) {
        return;
    }
    
    // 打开日志文件（追加模式）
    std::ofstream logFile;
    logFile.open(m_flowTableLogFilename.c_str(), std::ios::app | std::ios::out);
    
    if (!logFile.is_open()) {
        NS_LOG_ERROR("Failed to open flow table log file: " << m_flowTableLogFilename);
        return;
    }
    
    // 记录当前时间
    logFile << "\n==== Flow Table for Switch " << GetId() << " at " 
            << Simulator::Now().GetSeconds() << "s ====\n";
    logFile << "Total Active Flows: " << m_flowTable.size() << std::endl;
    
    // 记录流的详细信息
    uint64_t totalBytes = 0;
    uint64_t totalPackets = 0;
    uint64_t activeQueueBytes = 0;
    uint64_t activeQueuePackets = 0;
    
    for (const auto& entry : m_flowTable) {
        const FlowKey& key = entry.first;
        const FlowQueueStats& stats = entry.second;
        
        // 记录流的基本信息
        logFile << "Flow: " 
                << Ipv4Address(key.sip) << ":" << key.sport << " -> " 
                << Ipv4Address(key.dip) << ":" << key.dport 
                << " [Protocol: " << (int)key.protocol << "]\n";
                
        // 记录流的队列状态
        logFile << "  Queue: " << stats.queueBytes << " bytes, " 
                << stats.queuePackets << " packets\n";
        logFile << "  Total: " << stats.totalBytes << " bytes, " 
                << stats.totalPackets << " packets\n";
        logFile << "  Last Updated: " << stats.lastUpdateTime.GetSeconds() << "s\n";
        
        // 更新统计数据
        totalBytes += stats.totalBytes;
        totalPackets += stats.totalPackets;
        activeQueueBytes += stats.queueBytes;
        activeQueuePackets += stats.queuePackets;
    }
    
    // 记录汇总信息
    logFile << "\nSummary Statistics:\n";
    logFile << "  Active Queue: " << activeQueueBytes << " bytes, " 
            << activeQueuePackets << " packets\n";
    logFile << "  Cumulative: " << totalBytes << " bytes, " 
            << totalPackets << " packets\n";
    logFile << "  Avg Bytes Per Flow: " << (m_flowTable.empty() ? 0 : totalBytes / m_flowTable.size()) << "\n";
    logFile << "===========================\n";
    
    logFile.close();
    
    // 安排下一次记录
    m_flowTableLogEvent = Simulator::Schedule(MicroSeconds(m_flowTableLogInterval), 
                                              &SwitchNode::LogFlowTablePeriodically, this);
}
/**	Flow Table Logging**/


/** Swtich Node DCQCN **/
#define PRINT_LOG 0
/******************************
 * 
 * Mellanox's version of DCQCN
 *****************************/

void SwitchNode::UpdateAlphaMlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;

        auto mlx = dcqMap.find(key);

        #if PRINT_LOG
        //std::cout << Simulator::Now() << " alpha update:" << m_node->GetId() << ' ' << q->mlx.m_alpha << ' ' << (int)q->mlx.m_alpha_cnp_arrived << '\n';
        //printf("%lu alpha update: %08x %08x %u %u %.6lf->", Simulator::Now().GetTimeStep(), q->sip.Get(), q->dip.Get(), q->sport, q->dport, q->mlx.m_alpha);
        #endif

        // std::cout << "****" << std::endl;

        if (mlx->second.m_alpha_cnp_arrived){
                mlx->second.m_alpha = (1 - m_g)*mlx->second.m_alpha + m_g;      //binary feedback
        }else {
                mlx->second.m_alpha = (1 - m_g)*mlx->second.m_alpha;    //binary feedback
        }
        #if PRINT_LOG
        //printf("%.6lf\n", q->mlx.m_alpha);
        #endif
        mlx->second.m_alpha_cnp_arrived = false; // clear the CNP_arrived bit

//      dcqMap[key] = mlx;

        ScheduleUpdateAlphaMlx(key);
}

void SwitchNode::ScheduleUpdateAlphaMlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

        mlx->second.m_eventUpdateAlpha = Simulator::Schedule(MicroSeconds(m_alpha_resume_interval), &SwitchNode::UpdateAlphaMlx, this, key);

//      dcqMap[key] = mlx;
}

void SwitchNode::cnp_received_mlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

//      std::cout << mlx.m_first_cnp << std::endl;

        mlx->second.m_alpha_cnp_arrived = true; // set CNP_arrived bit for alpha update
        mlx->second.m_decrease_cnp_arrived = true; // set CNP_arrived bit for rate decrease
        if (mlx->second.m_first_cnp){
                // init alpha
                mlx->second.m_alpha = 1;
                mlx->second.m_alpha_cnp_arrived = false;
                // schedule alpha update
                ScheduleUpdateAlphaMlx(key);
                // schedule rate decrease
                ScheduleDecreaseRateMlx(key, 1); // add 1 ns to make sure rate decrease is after alpha update
                // set rate on first CNP
                mlx->second.m_targetRate = mlx->second.m_rate = m_rateOnFirstCNP * mlx->second.m_rate;
                mlx->second.m_first_cnp = false;
        }

//      dcqMap[key] = mlx;
}

void SwitchNode::CheckRateDecreaseMlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

//      std::cout << dcqMap[key].m_decrease_cnp_arrived  << std::endl;

        ScheduleDecreaseRateMlx(key, 0);

        if (mlx->second.m_decrease_cnp_arrived){
                bool clamp = true;
                if (!m_EcnClampTgtRate){
                        if (mlx->second.m_rpTimeStage == 0)
                                clamp = false;
                }

                if (clamp){
//                      std::cout << "XJXJXJXJXJX" << std::endl;
                        mlx->second.m_targetRate = mlx->second.m_rate;
                }

                mlx->second.m_rate = std::max(m_minRate, mlx->second.m_rate * (1 - mlx->second.m_alpha / 2));
                // std::cout << mlx.m_rate  << std::endl;

                // reset rate increase related things
                mlx->second.m_rpTimeStage = 0;
                mlx->second.m_decrease_cnp_arrived = false;
                Simulator::Cancel(mlx->second.m_rpTimer);

//              std::cout << "XJXJXJXJXJX" << std::endl;
                mlx->second.m_rpTimer = Simulator::Schedule(MicroSeconds(m_rpgTimeReset), &SwitchNode::RateIncEventTimerMlx, this, key);
        }

//      dcqMap[key] = mlx;

//      std::cout << dcqMap[key].m_decrease_cnp_arrived  << std::endl;
}

void SwitchNode::ScheduleDecreaseRateMlx(FlowKey key, uint32_t delta){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

        mlx->second.m_eventDecreaseRate = Simulator::Schedule(MicroSeconds(m_rateDecreaseInterval) + NanoSeconds(delta), &SwitchNode::CheckRateDecreaseMlx, this, key);

//      dcqMap[key] = mlx;
}

void SwitchNode::RateIncEventTimerMlx(FlowKey key){
//      std::cout << "XJXJXJXJXJX" << std::endl;

        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

        mlx->second.m_rpTimer = Simulator::Schedule(MicroSeconds(m_rpgTimeReset), &SwitchNode::RateIncEventTimerMlx, this, key);
//      std::cout << key  << std::endl;
        RateIncEventMlx(key);
        mlx->second.m_rpTimeStage++;

//      dcqMap[key] = mlx;
}

void SwitchNode::RateIncEventMlx(FlowKey key){
//      std::cout << key  << std::endl;

        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

//      std::cout << mlx.m_rpTimeStage  << std::endl;

        // check which increase phase: fast recovery, active increase, hyper increase
        if (mlx->second.m_rpTimeStage < m_rpgThreshold){ // fast recovery
                FastRecoveryMlx(key);
        }else if (mlx->second.m_rpTimeStage == m_rpgThreshold){ // active increase
//              std::cout << mlx.m_targetRate  << std::endl;
                ActiveIncreaseMlx(key);
        }else { // hyper increase
                HyperIncreaseMlx(key);
        }

        //dcqMap[key] = mlx;
}

void SwitchNode::FastRecoveryMlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

        mlx->second.m_rate = (mlx->second.m_rate / 2) + (mlx->second.m_targetRate / 2);
        //dcqMap[key] = mlx;
}

void SwitchNode::ActiveIncreaseMlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

        /* get NIC
        uint32_t nic_idx = GetNicIdxOfQp(q);
        Ptr<QbbNetDevice> dev = m_nic[nic_idx].dev;
        */
        // increate rate
        mlx->second.m_targetRate += m_rai;
//      std::cout << mlx.m_targetRate  << std::endl;

        if (mlx->second.m_targetRate > lineRate)
                mlx->second.m_targetRate = lineRate;

        mlx->second.m_rate = (mlx->second.m_rate / 2) + (mlx->second.m_targetRate / 2);

//      dcqMap[key] = mlx;
}

void SwitchNode::HyperIncreaseMlx(FlowKey key){
        if(dcqMap.find(key) == dcqMap.end()) return;
        auto mlx = dcqMap.find(key);

        /* get NIC
        uint32_t nic_idx = GetNicIdxOfQp(q);
        Ptr<QbbNetDevice> dev = m_nic[nic_idx].dev;
        */

        // increate rate

        mlx->second.m_targetRate += m_rhai;
        if (mlx->second.m_targetRate > lineRate)
                mlx->second.m_targetRate = lineRate;

        mlx->second.m_rate = (mlx->second.m_rate / 2) + (mlx->second.m_targetRate / 2);

//      dcqMap[key] = mlx;
}
/** Swtich Node DCQCN **/

} /* namespace ns3 */
