#ifndef SWITCH_NODE_H
#define SWITCH_NODE_H

#include <unordered_map>
#include <ns3/node.h>
#include "qbb-net-device.h"
#include "switch-mmu.h"
#include "pint.h"
#include <ns3/event-id.h>

/** DCI NUMBER **/
# define DCI_SWITCH_0 0
# define DCI_SWITCH_1 37
/** DCI NUMBER **/

namespace ns3 {

/** DCQCN **/
struct Mlx {
        DataRate m_rate;

        DataRate m_targetRate;  //< Target rate
        EventId m_eventUpdateAlpha;
        double m_alpha;
        bool m_alpha_cnp_arrived; // indicate if CNP arrived in the last slot
        bool m_first_cnp; // indicate if the current CNP is the first CNP
        EventId m_eventDecreaseRate;
        bool m_decrease_cnp_arrived; // indicate if CNP arrived in the last slot
        uint32_t m_rpTimeStage;
        EventId m_rpTimer;
};
/** DCQCN **/

/*
struct FlowKey {
        uint32_t sip;          // 源IP
        uint32_t dip;          // 目的IP
        uint16_t sport;        // 源端口
        uint16_t dport;        // 目的端口
        uint8_t protocol;      // 协议类型

        // 用于哈希表的相等比较
        bool operator==(const FlowKey &other) const {
       		 return sip == other.sip && dip == other.dip &&
                        sport == other.sport && dport == other.dport &&
                        protocol == other.protocol;
	}
};

        // 自定义哈希函数
        struct FlowKeyHash {
        std::size_t operator()(const FlowKey& k) const {
                return std::hash<uint32_t>()(k.sip) ^
                   std::hash<uint32_t>()(k.dip) ^
                   std::hash<uint16_t>()(k.sport) ^
                   std::hash<uint16_t>()(k.dport) ^
                  std::hash<uint8_t>()(k.protocol);
        }
        };

        // 每个流的队列统计信息
        struct FlowQueueStats {
                uint32_t queueBytes;       // 当前队列中的字节数
                uint32_t queuePackets;     // 当前队列中的数据包数
                uint32_t maxQueueBytes;    // 历史最大队列字节数
        uint32_t maxQueuePackets;  // 历史最大队列包数
                uint32_t outDev;           // 出口设备索引
                uint32_t qIndex;           // 队列索引/优先级
                uint64_t totalBytes;       // 累计字节数
                uint64_t totalPackets;     // 累计数据包数
                Time lastUpdateTime;       // 最后更新时间
        };
        // 流表类型定义
        typedef std::unordered_map<FlowKey, FlowQueueStats, FlowKeyHash> FlowTable;

*/

class Packet;

class SwitchNode : public Node{
	static const uint32_t pCnt = 257;	// Number of ports used
	static const uint32_t qCnt = 8;	// Number of queues/priorities used
	uint32_t m_ecmpSeed;
	std::unordered_map<uint32_t, std::vector<int> > m_rtTable; // map from ip address (u32) to possible ECMP port (index of dev)

	// monitor of PFC
	uint32_t m_bytes[pCnt][pCnt][qCnt]; // m_bytes[inDev][outDev][qidx] is the bytes from inDev enqueued for outDev at qidx
	
	uint64_t m_txBytes[pCnt]; // counter of tx bytes

	uint32_t m_lastPktSize[pCnt];
	uint64_t m_lastPktTs[pCnt]; // ns
	double m_u[pCnt];

protected:
	bool m_ecnEnabled;
	uint32_t m_ccMode;
	uint64_t m_maxRtt;

	uint32_t m_ackHighPrio; // set high priority for ACK/NACK

private:
	int GetOutDev(Ptr<const Packet>, CustomHeader &ch);
	void SendToDev(Ptr<Packet>p, CustomHeader &ch);
	static uint32_t EcmpHash(const uint8_t* key, size_t len, uint32_t seed);
	void CheckAndSendPfc(uint32_t inDev, uint32_t qIndex);
	void CheckAndSendResume(uint32_t inDev, uint32_t qIndex);
public:
	Ptr<SwitchMmu> m_mmu;
	/** test **/
	uint64_t lastBytes;

	static TypeId GetTypeId (void);
	SwitchNode();
	void SetEcmpSeed(uint32_t seed);
	void AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx);
	void ClearTable();
	bool SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch);
	void SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p);
	
	/** test **/
	EventId m_rpTimer;
	double t1;
	void TimerEvent(uint32_t inDev, uint32_t qIndex);

	// for approximate calc in PINT
	int logres_shift(int b, int l);
	int log2apprx(int x, int b, int m, int l); // given x of at most b bits, use most significant m bits of x, calc the result in l bits

	/** BICC **/
	bool isSendDCI;
	bool isRecvDCI;
	uint64_t maxBW;
	void sendCNPByDCI(Ptr<Packet> p, uint32_t ifIndex);
	
	std::unordered_map<uint32_t, uint64_t> forward_table;
	std::unordered_map<uint32_t, uint64_t> recv_table;
	uint64_t iBDP;
	// std::unordered_map<uint64_t, uint64_t> through_table;
	// std::unordered_map<FlowKey, uint64_t, FlowKeyHash> through_table2;

	EventId rpTimer;
	uint64_t counter;
	double TimeReset;
	void CalcEvent();
	/** BiCC **/

	/**	Flow Table **/

	struct FlowKey {
		uint32_t sip;          // 源IP
		uint32_t dip;          // 目的IP
		uint16_t sport;        // 源端口
		uint16_t dport;        // 目的端口
		uint8_t protocol;      // 协议类型
		
		// 用于哈希表的相等比较
		bool operator==(const FlowKey &other) const {
			return sip == other.sip && dip == other.dip && 
				   sport == other.sport && dport == other.dport && 
				   protocol == other.protocol;
		}
	};

	// 自定义哈希函数
	struct FlowKeyHash {
    	std::size_t operator()(const FlowKey& k) const {
        	return std::hash<uint32_t>()(k.sip) ^ 
            	   std::hash<uint32_t>()(k.dip) ^ 
            	   std::hash<uint16_t>()(k.sport) ^ 
               	   std::hash<uint16_t>()(k.dport) ^
                  std::hash<uint8_t>()(k.protocol);
    	}
	};

	// 每个流的队列统计信息
	struct FlowQueueStats {
		uint32_t queueBytes;       // 当前队列中的字节数
		uint32_t queuePackets;     // 当前队列中的数据包数
		uint32_t maxQueueBytes;    // 历史最大队列字节数 
    	uint32_t maxQueuePackets;  // 历史最大队列包数 
		uint32_t outDev;           // 出口设备索引
		uint32_t qIndex;           // 队列索引/优先级
		uint64_t totalBytes;       // 累计字节数
		uint64_t totalPackets;     // 累计数据包数
		Time lastUpdateTime;       // 最后更新时间
		Time lastCnpTime;  		   // 该流最后 CNP 的发送时间
		/** rate monitor **/
		uint64_t currentRate; // 当前速率
		/** rate monitor **/
	};
	// 流表类型定义
	typedef std::unordered_map<FlowKey, FlowQueueStats, FlowKeyHash> FlowTable;

	// 流表，存储每个流的统计信息
	FlowTable m_flowTable; 

	// 从包和头部提取FlowKey

    	FlowKey ExtractFlowKey(CustomHeader &ch);

	void PrintFlowTable(bool detailed = false);

	// 流表日志相关
	bool m_flowTableLoggingEnabled;

	uint32_t m_targetFlowTableSwitchId;
   	std::string m_flowTableLogFilename;
    	double m_flowTableLogInterval;
    	EventId m_flowTableLogEvent;

	void LogFlowTablePeriodically();
   	void StartFlowTableLogging(uint32_t targetSwitchId, double intervalMs, std::string filename);
	/**	Flow Table **/

	/** Control Message **/
	uint32_t queueThreshold; // 队列阈值
	bool DCI_CM_test;

	void SendControlMessage(Ptr<Packet> p, uint32_t ifIndex);
	/** Control Message **/

	/** SRC-DCI CNP GENERATION **/
	void SrcDCICNPGen(Ptr<Packet> p, uint32_t ifIndex);
	/** SRC-DCI CNP GENERATION **/


	/** DCI_ALG **/
	bool m_dciAlgEnabled;
	/** DCI_ALG **/

	std::unordered_map<FlowKey, uint64_t, FlowKeyHash> through_table2;

    	
        /** DCQCN And RateLimiter **/
        /******************************
         * Mellanox's version of DCQCN
         *****************************/
        double m_g; //feedback weight
        double m_rateOnFirstCNP; // the fraction of line rate to set on first CNP
        bool m_EcnClampTgtRate;
        double m_rpgTimeReset;
        double m_rateDecreaseInterval;
        uint32_t m_rpgThreshold;
        double m_alpha_resume_interval;
        DataRate m_rai;         //< Rate of additive increase
        DataRate m_rhai;                //< Rate of hyper-additive increase

        // the Mellanox's version of alpha update:
        // every fixed time slot, update alpha.
        void UpdateAlphaMlx(FlowKey key);
        void ScheduleUpdateAlphaMlx(FlowKey key);

        // Mellanox's version of CNP receive
        void cnp_received_mlx(FlowKey key);

        // Mellanox's version of rate decrease
        // It checks every m_rateDecreaseInterval if CNP arrived (m_decrease_cnp_arrived).
        // If so, decrease rate, and reset all rate increase related things
        void CheckRateDecreaseMlx(FlowKey key);
        void ScheduleDecreaseRateMlx(FlowKey key, uint32_t delta);

        // Mellanox's version of rate increase
        DataRate m_minRate;             //< Min sending rate
        DataRate lineRate;
        uint32_t slidingWin;
        uint32_t initWin;
        // std::unordered_map<uint64_t, Mlx> dcqMap;
	std::unordered_map<FlowKey, Mlx, FlowKeyHash> dcqMap;

        void RateIncEventTimerMlx(FlowKey key);
        void RateIncEventMlx(FlowKey key);
        void FastRecoveryMlx(FlowKey key);
        void ActiveIncreaseMlx(FlowKey key);
        void HyperIncreaseMlx(FlowKey key);

        /******************************
         * Token Bucket Rate Limiting
         *****************************/
//        std::unordered_map<uint64_t, uint64_t> tokenBuckets;
   	std::unordered_map<FlowKey, uint64_t, FlowKeyHash> tokenBuckets;

   	/** DCQCN And RateLimiter **/

	/** **/
	 struct QueuedPacketInfo {
		Ptr<Packet> packet;
		int outDevIdx;
		uint32_t qIndex; 
	};

	std::unordered_map<FlowKey,std::queue<QueuedPacketInfo>,FlowKeyHash> m_perFlowQueues;

};

} /* namespace ns3 */

#endif /* SWITCH_NODE_H */
