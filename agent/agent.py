from __future__ import absolute_import, division, print_function

from chunk import Chunk
import io
import json
import struct
from confluent_kafka import Producer, Consumer, KafkaError

import agent.event as event


# Chunk header: type and size.
# TODO(jerry): Move CHUNK_HEADER to a shared module.
CHUNK_HEADER = struct.Struct('>4s I')
JSON_CHUNK_TYPE = b'JSON'  # a JSON chunk contains a JSON dictionary message
BLOB_CHUNK_TYPE = b'BLOB'  # a BLOB chunk contains a bytes

BLOBS = 'blobs' # message dict key for a sequence of binary large objects
				# (bytes), unsuitable for the JSON part of the payload, or to
				# print to a log, or to put in a command-line arg.


def delivery_report(err, msg):
	"""
	This is a utility method passed to the Kafka Producer to handle the delivery
	of messages sent using `send(topic, message)`. 
	"""
	if err is not None:
		print('message delivery failed: {}'.format(msg))
		print('failed message: {}'.format(err))

class Agent(object):

	"""
	Agent: a process with the ability to send and receive messages through Kafka.

	This is a base class that handles all of the details of connecting to Kafka and the 
	initialization and configuration of a Producer (to send messages) and Consumer
	(to receive messages). 

	To subclass Agent:

	* First, a call to super must be made in the `__init__` method of the overriding class.
	  `agent_id` is any unique string and identifies this agent with other agents in the system.
	  `kafka_config` is a configuration dictionary containing information about the kafka host,
	  topics to subscribe to and additional information. Its details are described in the
	  docstring for `__init__`.

	* Then, you can override `initialize()` to send any messages before the consumer
	  loop begins.

	* The main method to override is `receive(topic, message)`, which will be called
	  each time a message is received on any of the topics originally subscribed to.

	* Finally, on shutdown there is a call to `finalize()` where you can put anything
	  that needs to be cleaned up or any final messages to send before the process is
	  terminated.

	Additionally, a `send(topic, message)` method is provided that does not need to be
	overridden but can be used directly, which sends a message to the specified topic
	using the configured Kafka Producer.
	"""

	def __init__(self, agent_id, agent_type, agent_config):
		"""
		Initialize the Agent object with a unique id and kafka configuration.

		Args:
		    agent_id (str): A unique identifier which distinguishes this agent from
		        the rest of the agents in the system.

		    agent_type (str): A string indicating which type of agent this is. This is helpful
		        for spawning new agents of the same type from the agent shepherd.

		    agent_config (dict): A dictionary containing any configuration information the agent
		        might need. Subclasses use this extensively, but the only key that needs to be
		        present for the base agent class is:

		        * kafka_config (dict): A dictionary containing all of the kafka configuration 
                    information, including:

		            * host (str): The address of the host machine running Kakfa.
		            * subscribe (array[str]): A list of topics for the Consumer to subscribe. 
		                to. If this is an empty array then no Consumer will be initialized.
		            * topics (dict[str, str]): A mapping from topic roles to specific topics.
		                These topics provide a way to send messages to the various agents in the
		                system.
		"""

		self.agent_id = agent_id
		self.agent_type = agent_type
		self.agent_config = agent_config
		self.kafka_config = agent_config['kafka_config']
		self.topics = self.kafka_config['topics']

		self.producer = Producer({
			'bootstrap.servers': self.kafka_config['host']})

		self.running = False
		self.initialized = False
		self.consumer = None
		if self.kafka_config['subscribe']:
			self.consumer = Consumer({
				'bootstrap.servers': self.kafka_config['host'],
				'enable.auto.commit': True,
				'group.id': 'simulation-' + str(agent_id),
				'default.topic.config': {
					'auto.offset.reset': 'latest'}})

	def preinitialize(self):
		"""
		Initialize the Agent in the system.

		This method is called after the Producer and Consumer have been initialized.
		It is a good place to send any initialization messages to other agents before
		falling into the polling cycle.
		"""

		pass

	def start(self):
		"""
		Start the polling loop if we have a consumer, otherwise just call `initialize` directly
		"""

		if self.consumer:
			topics = self.kafka_config['subscribe'] + [self.topics['agent_receive']]
			self.consumer.subscribe(topics)

			self.poll()
		else:
			self.preinitialize()
			self.initialized = True

	def poll(self):
		"""
		Enter the main consumer polling loop.

		Once poll is called, the thread will be claimed and any interaction with the 
		system from this point on will be mediated through message passing. This is called
		at the end of the base class's `__init__(agent_id, kafka_config)` method and does not
		need to be called manually by the subclass.
		"""

		self.running = True
		while self.running:
			raw = self.consumer.poll(timeout=1.0)  # timeout (in seconds) so ^C works

			# calling initialize() once consumer is established so as not to miss
			# immediate responses to initialization sends. If `poll` is not called before an
			# initialization message is sent then an immediate response could be missed.
			if not self.initialized:
				self.preinitialize()
				self.initialized = True

			if raw is None:
				continue
			if raw.error():
				if raw.error().code() == KafkaError._PARTITION_EOF:
					continue
				else:
					print('Error in kafka consumer:', raw.error())

					self.running = False

			else:
				# `raw.value()` is implemented in C with a docstring that
				# suggests it needs a `payload` argument. Suppress the warning.
				# noinspection PyArgumentList
				message = self.decode_payload(raw.value())
				if not message:
					continue

				if message['event'] == event.GLOBAL_SHUTDOWN:
					self.shutdown()
				else:
					self.receive(raw.topic(), message)

	def send(self, topic, message, print_send=True):
		"""
		Send a Kafka message on the given topic. The payload is a transmitted
		as a JSON dictionary chunk and optional BLOB chunks.

		Args:
			topic (str): The Kafka topic to send the message on.

			message (dict): A dictionary containing the message to send. This dictionary
				needs to be JSON serializable, so it must contain only basic types like `str`,
				`int`, `float`, `list`, `tuple`, `array`, and `dict`. Any functions or objects
				present will throw errors.

				SPECIAL CASE: `message[BLOBS]`, if present, is a sequence of BLOBs (bytes).

		    print_send (bool): Whether or not to print the message that is sent.
		"""
		payload = self.encode_payload(message)

		if print_send:
			self.print_message(topic, message, False)

		self.producer.produce(
			topic,
			payload,
			callback=delivery_report)

		self.producer.flush(timeout=1.0)

	def print_message(self, topic, message, incoming=True):
		"""Print the incoming or outgoing message to the console/log, redacting
		any BLOBs.
		"""
		num_blobs = 0
		redacted_message = dict(message)
		if BLOBS in message:
			num_blobs = len(message[BLOBS])
			del redacted_message[BLOBS]

		print('{} {} {}'       # <-- topic event
			  ' [{} {}]:'      # [agent_type agent_id]
			  ' {}{}'.format(  # {message dict} + 2 BLOBs
			'-->' if incoming else '<--',
			topic,
			message.get('event', 'generic'),

			self.agent_type,
			self.agent_id,

			redacted_message,
			' + {} BLOBs'.format(num_blobs) if num_blobs else '',
			))

	def encode_payload(self, message):
		"""Encode a `message` dictionary as a Kafka message payload chunk
		stream. The first chunk holds the JSON message. Following chunks hold
		BLOBs extracted from `message[BLOBS]`.
		"""
		message = dict(message)
		blobs = message.pop(BLOBS, [])

		# json.dumps(m, ensure_ascii=False) returns a str or unicode string, depending on
		# content (always a unicode string in Python 3) w/o \u escapes. Encode that into
		# UTF-8 bytes. print() can decode UTF-8 bytes but not \u escapes.
		encoded = json.dumps(message, ensure_ascii=False).encode('utf-8')

		payload = [
			CHUNK_HEADER.pack(JSON_CHUNK_TYPE, len(encoded)),
			encoded]
		for blob in blobs:
			payload.append(CHUNK_HEADER.pack(BLOB_CHUNK_TYPE, len(blob)))
			payload.append(blob)

		return b''.join(payload)

	def decode_payload(self, payload):
		"""Decode a Kafka message payload chunk stream."""
		bytestream = io.BytesIO(payload)

		chunk = Chunk(bytestream, align=False)
		if chunk.getname() != JSON_CHUNK_TYPE:
			print('Error: Each Agent message must start with a JSON chunk')
			return {}

		json_message = chunk.read().decode('utf-8')
		message = json.loads(json_message)
		chunk.close()

		blobs = []
		while True:
			try:
				chunk = Chunk(bytestream, align=False)
			except EOFError:
				break

			if chunk.getname() == BLOB_CHUNK_TYPE:
				blobs.append(chunk.read())

			chunk.close()

		message[BLOBS] = blobs
		return message

	def receive(self, topic, message):
		"""
		Receive a message from the given topic.

		This is the main method that needs to be overridden by a subclass to provide the 
		behavior for the agent. This method is invoked each time the agent receives a message
		on one of the topics originally subscribed to, which was supplied in the 
		kafka_config['subscribe'] array.

		By convention, each message contains an `event` key which can be switched on in this
		method to trigger a specific response.

		`message[BLOBS]` is a list (often empty) of bytes BLOBs.
		Don't print the BLOBs or send them on a command line. Call
		`self.print_message(topic, message)` to print/log the message.

		Args:
		    topic (str): The Kafka topic this message was received on.
		    message (dict): A deserialized dictionary containing at the least an `event` key
		        describing the event that originally triggered this message.
		"""

		pass

	def shutdown(self):
		"""
		Shutdown the Kafka Producer and Consumer.

		Breaks out of the polling receive loop and calls `finalize`, then flushes the
		producer, commits the consumer offsets and exits.
		"""

		self.running = False
		self.finalize()

		self.producer.flush(timeout=1.0)

		if self.consumer:
			self.consumer.commit()
			self.consumer.close()

	def finalize(self):
		"""
		Send any final messages and do any clean up before exiting.

		This method can be overridden by the subclass to send any final messages to other
		agents in the system or release system resources before the agent is shut down.
		"""

		pass
