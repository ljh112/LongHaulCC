/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2009 New York University
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
 * Author: Adrian S.-W. Tam <adrian.sw.tam@gmail.com>
 */

#include <stdint.h>
#include <iostream>
#include "ns3/buffer.h"
#include "ns3/address-utils.h"
#include "ns3/log.h"
#include "bicc-header.h"

NS_LOG_COMPONENT_DEFINE ("BICCHeader");

namespace ns3 {

NS_OBJECT_ENSURE_REGISTERED (BICCHeader);

BICCHeader::BICCHeader ()
{ 
  dport=0;
}


BICCHeader::~BICCHeader ()
{}

TypeId
BICCHeader::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::BICCHeader")
    .SetParent<Header> ()
    .AddConstructor<BICCHeader> ()
    ;
  return tid;
}

TypeId
BICCHeader::GetInstanceTypeId (void) const
{
  return GetTypeId ();
}

void BICCHeader::Print (std::ostream &os) const
{
  return;
}
uint32_t BICCHeader::GetSerializedSize (void)  const
{
  return sizeof(qIndex) + sizeof(dport);
}
void BICCHeader::Serialize (Buffer::Iterator start)  const
{
  start.WriteU16(qIndex);	
  start.WriteU16(dport);
}

uint32_t BICCHeader::Deserialize (Buffer::Iterator start)
{
  qIndex = start.ReadLsbtohU16();
  dport=start.ReadLsbtohU16();
  return GetSerializedSize ();
}


}; // namespace ns3
