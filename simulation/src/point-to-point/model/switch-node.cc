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
                        "The capacity of bandwitch on DCI switch for BiCC (Gbps)",
                        DoubleValue(200),
                        MakeDoubleAccessor(&SwitchNode::maxBW),
                        MakeDoubleChecker<uint64_t>())
			
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

	// BICC
//	isSendDCI = false;
//	isRecvDCI = false;
	// BICC
	
	forward_table.clear();
	recv_table.clear();

	through_table.clear();

	/** BICC **/
	rpTimer = Simulator::Schedule(MicroSeconds(TimeReset), &SwitchNode::CalcEvent, this);	
	counter = 0;	
	/** BICC **/

	/** Flow Table Logging **/
	/** 初始化 **/
	m_flowTableLoggingEnabled = false;
	m_targetFlowTableSwitchId = 0;
	m_flowTableLogInterval = 1000000.0; // 默认1秒
	m_flowTableLogFilename = "switch_flow_table.log";
	/** Flow Table Logging **/
}


void SwitchNode::CalcEvent(){
	if(37 == m_mmu->node_id){
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

                if(through_table.size() > 0){
//                      uint64_t val = std::min((uint64_t)(counter * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024), maxBW);
//                      uint64_t val = (uint64_t)(counter * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024);
//                      uint64_t val = counter;
//                      std::cout << clock() << " " << val << std::endl;
//                      printf("%u %u\n", clock(), val);

                        // through_table[ch.sip] += p->GetSize();
			uint64_t totalCnt = 0;
                        for(auto kv : through_table){
                                uint64_t key = kv.first;
                                uint64_t cnt = kv.second;
                                uint64_t val = std::min((uint64_t)(cnt * 1.0 * 8  / TimeReset * 1e6 /1024 / 1024 / 1024), maxBW);
				totalCnt += val;
//                                printf("%llu %llu %llu\n", clock(), key, val);
//  				std::cout << through_table.size() << std::endl;
  			}
			
			printf("%llu %llu\n", clock(), std::min(totalCnt,maxBW));
//			std::cout << through_table.size() << std::endl;
//                        through_table.clear();
                }


//             counter = 0;
		through_table.clear();


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

void SwitchNode::SendToDev(Ptr<Packet>p, CustomHeader &ch){
	/**	Flow Table**/
	/** 包到达处对对应流表项更新 **/
	if (ch.l3Prot == 0x11){ // 针对UDP流量进行信息统计
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

	int idx = GetOutDev(p, ch);
//	std::cout << "*****" << std::endl;
//	counter += p->GetSize();
	
//	std::cout << ch.sip << std::endl;
//	through_table[ch.sip] += p->GetSize();
//	std::cout << ch.sip  << std::endl;
        if (ch.l3Prot == 0x11 || ch.l3Prot == 0x6){
//              std::cout << ch.l3Prot << " "  << ch.sip << std::endl;
                // through_table[ch.sip] += p->GetSize();
		uint64_t flowId = ((uint64_t)ch.sip << 32) | ((uint64_t)ch.udp.pg << 16) | (uint64_t)ch.udp.sport;
		// std::cout << key << std::endl;
		through_table[flowId] += p->GetSize();
        }
	else if(ch.l3Prot == 0xFC || ch.l3Prot == 0xFD){ // || ch.l3Prot == 0xFF || ch.l3Prot == 0xFE
//              std::cout << ch.l3Prot << " "  << ch.dip << std::endl;
                // through_table[ch.dip] += p->GetSize();
                uint64_t flowId = ((uint64_t)ch.dip << 32) | ((uint64_t)ch.ack.pg << 16) | (uint64_t)ch.ack.dport;
                // std::cout << key << std::endl;
                through_table[flowId] += p->GetSize();
        }



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
			CheckAndSendPfc(inDev, qIndex);
		}

		/* BICC */
		uint8_t ecnbits = ch.GetIpv4EcnBits();
                bool egressCongested = m_ecnEnabled && m_mmu->ShouldSendCN(idx, qIndex);
                if(9 == m_ccMode && 0 == m_mmu->node_id && ch.l3Prot == 0x11 && (egressCongested || ecnbits)){
			sendCNPByDCI(p, idx);

			PppHeader ppp;
                        Ipv4Header h;
                        p->RemoveHeader(ppp);
                        p->RemoveHeader(h);
                        h.SetEcn((Ipv4Header::EcnType)0x00);
                        p->AddHeader(h);
                        p->AddHeader(ppp);
		}
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

		// counter += p->GetSize();		
//		m_devices[idx]->SwitchSend(qIndex, p, ch);
		
		/** BiCC third loop **/
		Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(m_devices[idx]);
		uint32_t fip = ch.dip;

		if(9 == m_ccMode && 37 == m_mmu->node_id && ch.l3Prot == 0x11){
			if(forward_table.find(fip)!=forward_table.end() 
					&& recv_table.find(fip)!=recv_table.end()
				      	&& forward_table[fip]>iBDP){
//				 std::cout << recv_table[fip] << " ******1******" << std::endl;
				return;
			}

			forward_table[fip] += p->GetSize();
//			std::cout << bi_table[fip] << " ******1******" << std::endl;
		}	
	
		if(9 == m_ccMode && 37 == m_mmu->node_id && ch.l3Prot == 0xFC){
//			std::cout << ch.udp.ih.recvBiCCBytes << " ******1******" << std::endl;
//			uint32_t aSip = ch.sip;
			recv_table[ch.sip] = ch.ack.ih.recvBiCCBytes;			
//			if(bi_table[aSip] < ch.udp.ih.recvBiCCBytes)
//			std::cout << forward_table[aSip]  << " ******2******" << recv_table[aSip]  << std::endl;
//			if(bi_table.find(fip) != bi_table.end()){
//				bi_table[fip] -= ch.udp.ih.recvBiCCBytes;
//				std::cout << bi_table[fip] << std::endl;
//				
//				if(bi_table[fip] <= 0){
//					bi_table.erase(fip);
//				}
//			}
//			ch.udp.ih.recvBiCCBytes = 1;
//			std::cout << ch.udp.ih.isLongLoop << " ******2******" << std::endl;	
		}

		/** BiCC third loop **/	
		
		m_devices[idx]->SwitchSend(qIndex, p, ch);

	//	std::cout << dev->GetDataRate().GetBitRate() << std::endl;
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
//	std::cout << "*****" << std::endl;

	FlowIdTag t;
	p->PeekPacketTag(t);
	counter += p->GetSize(); // test
///        through_table[ch.sip] += p->GetSize();	
//	uint32_t sip = 0x0b000001 + ((m_id / 256) * 0x00010000) + ((m_id % 256) * 0x00000100);
	CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
        ch.getInt = 1; // parse INT header
        p->PeekHeader(ch);
//	through_table[ch.sip] += p->GetSize();
	
	if (ch.l3Prot == 0x11 || ch.l3Prot == 0x6){
//		std::cout << ch.l3Prot << " "  << ch.sip << std::endl;	
//		through_table[ch.sip] += p->GetSize();
		uint64_t flowId = ((uint64_t)ch.sip << 32) | ((uint64_t)ch.udp.pg << 16) | (uint64_t)ch.udp.sport;
                through_table[flowId] += p->GetSize();
	}
	else if(ch.l3Prot == 0xFC || ch.l3Prot == 0xFD){ // || ch.l3Prot == 0xFF || ch.l3Prot == 0xFE
//		std::cout << ch.l3Prot << " "  << ch.dip << std::endl;
//		through_table[ch.dip] += p->GetSize();
                uint64_t flowId = ((uint64_t)ch.dip << 32) | ((uint64_t)ch.ack.pg << 16) | (uint64_t)ch.ack.dport;
                through_table[flowId] += p->GetSize();
	}

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
			//std::cout << "22222222222" << std::endl;

			if (egressCongested){
				PppHeader ppp;
				Ipv4Header h;
				p->RemoveHeader(ppp);
				p->RemoveHeader(h);
				h.SetEcn((Ipv4Header::EcnType)0x03);
				p->AddHeader(h);
				p->AddHeader(ppp);

				/** BICC **/
//				std::cout << m_mmu->node_id << std::endl;
//				if(0 == m_mmu->node_id){
//					
//				}

				/** BICC **/
			}
		}
		//CheckAndSendPfc(inDev, qIndex);
                
		/** test **/
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

//			if(lastBytes>0 && lastBytes>m_mmu->GetSharedUsed(inDev, qIndex)){
//				std::cout << inDev << " " << qIndex << " "  << m_mmu->paused[inDev][qIndex] << " " << m_mmu->GetSharedUsed(inDev, qIndex) << " " << m_mmu->hdrm_bytes[inDev][qIndex] << std::endl;
                                Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
                                
				if(true == m_mmu->paused[inDev][qIndex]){
					device->SendPfc(qIndex, 1);	
				}

				m_mmu->SetResume(inDev, qIndex);
//			}
			
//			lastBytes = m_mmu->GetSharedUsed(inDev, qIndex);
		}
		*/	
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

} /* namespace ns3 */
