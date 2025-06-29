#ifndef CondCore_CondDB_IOVSchema_h
#define CondCore_CondDB_IOVSchema_h

#include "DbCore.h"
#include "IDbSchema.h"
//
#include <boost/date_time/posix_time/posix_time.hpp>

namespace cond {

  namespace persistency {

    conddb_table(TAG) {
      conddb_column(NAME, std::string);
      conddb_column(TIME_TYPE, cond::TimeType);
      conddb_column(OBJECT_TYPE, std::string);
      conddb_column(SYNCHRONIZATION, cond::SynchronizationType);
      conddb_column(END_OF_VALIDITY, cond::Time_t);
      conddb_column(DESCRIPTION, std::string);
      conddb_column(LAST_VALIDATED_TIME, cond::Time_t);
      conddb_column(INSERTION_TIME, boost::posix_time::ptime);
      conddb_column(MODIFICATION_TIME, boost::posix_time::ptime);
      conddb_column(PROTECTION_CODE, int);

      class Table : public ITagTable {
      public:
        explicit Table(coral::ISchema& schema);
        ~Table() override {}
        bool exists() override;
        void create() override;
        bool select(const std::string& name) override;
        bool select(const std::string& name,
                    cond::TimeType& timeType,
                    std::string& objectType,
                    cond::SynchronizationType& synchronizationType,
                    cond::Time_t& endOfValidity,
                    cond::Time_t& lastValidatedTime,
                    int& protectionCode) override;
        bool getMetadata(const std::string& name,
                         std::string& description,
                         boost::posix_time::ptime& insertionTime,
                         boost::posix_time::ptime& modificationTime) override;
        void insert(const std::string& name,
                    cond::TimeType timeType,
                    const std::string& objectType,
                    cond::SynchronizationType synchronizationType,
                    cond::Time_t endOfValidity,
                    const std::string& description,
                    cond::Time_t lastValidatedTime,
                    const boost::posix_time::ptime& insertionTime) override;
        void update(const std::string& name,
                    cond::SynchronizationType synchronizationType,
                    cond::Time_t& endOfValidity,
                    cond::Time_t lastValidatedTime,
                    const boost::posix_time::ptime& updateTime) override;
        void updateMetadata(const std::string& name,
                            const std::string& description,
                            const boost::posix_time::ptime& updateTime) override;
        void updateValidity(const std::string& name,
                            cond::Time_t lastValidatedTime,
                            const boost::posix_time::ptime& updateTime) override;
        void setProtectionCode(const std::string& name, int code) override;
        void unsetProtectionCode(const std::string& name, int code) override;

        bool isProtectable() { return m_isProtectable; }

      private:
        coral::ISchema& m_schema;
        bool m_isProtectable = false;
      };
    }

    conddb_table(PAYLOAD) {
      static constexpr unsigned int PAYLOAD_HASH_SIZE = 40;

      conddb_column(HASH, std::string, PAYLOAD_HASH_SIZE);
      conddb_column(OBJECT_TYPE, std::string);
      conddb_column(DATA, cond::Binary);
      conddb_column(STREAMER_INFO, cond::Binary);
      conddb_column(VERSION, std::string);
      conddb_column(INSERTION_TIME, boost::posix_time::ptime);

      class Table : public IPayloadTable {
      public:
        explicit Table(coral::ISchema& schema);
        ~Table() override {}
        bool exists() override;
        void create() override;
        bool select(const cond::Hash& payloadHash);
        bool select(const cond::Hash& payloadHash,
                    std::string& objectType,
                    cond::Binary& payloadData,
                    cond::Binary& streamerInfoData) override;
        bool getType(const cond::Hash& payloadHash, std::string& objectType) override;
        bool insert(const cond::Hash& payloadHash,
                    const std::string& objectType,
                    const cond::Binary& payloadData,
                    const cond::Binary& streamerInfoData,
                    const boost::posix_time::ptime& insertionTime);
        cond::Hash insertIfNew(const std::string& objectType,
                               const cond::Binary& payloadData,
                               const cond::Binary& streamerInfoData,
                               const boost::posix_time::ptime& insertionTime) override;

      private:
        coral::ISchema& m_schema;
      };
    }

    conddb_table(IOV) {
      conddb_column(TAG_NAME, std::string);
      conddb_column(SINCE, cond::Time_t);
      conddb_column(PAYLOAD_HASH, std::string, PAYLOAD::PAYLOAD_HASH_SIZE);
      conddb_column(INSERTION_TIME, boost::posix_time::ptime);
      static constexpr auto minSINCE_ =
          condcore_detail::addMin<SINCE::fullyQualifiedName().size()>(SINCE::fullyQualifiedName());
      static constexpr std::string_view minSINCE() { return std::string_view(minSINCE_.data()); }
      static constexpr auto maxSINCE_ =
          condcore_detail::addMax<SINCE::fullyQualifiedName().size()>(SINCE::fullyQualifiedName());
      static constexpr std::string_view maxSINCE() { return std::string_view(maxSINCE_.data()); }

      struct SINCE_GROUP {
        typedef cond::Time_t type;
        static constexpr size_t size = 0;
        static constexpr std::string_view tableName() { return SINCE::tableName(); }
        static constexpr std::string_view fullyQualifiedName() { return minSINCE(); }
        static std::string group(unsigned long long groupSize) {
          std::string sgroupSize = std::to_string(groupSize);
          return "CAST(" + std::string(SINCE::fullyQualifiedName()) + "/" + sgroupSize + " AS INT )*" + sgroupSize;
        }
      };

      struct SEQUENCE_SIZE {
        typedef unsigned int type;
        static constexpr size_t size = 0;
        static constexpr std::string_view tableName() { return SINCE::tableName(); }
        static constexpr std::string_view fullyQualifiedName() { return "COUNT(*)"; }
      };

      struct MIN_SINCE {
        typedef cond::Time_t type;
        static constexpr size_t size = 0;
        static constexpr std::string_view tableName() { return SINCE::tableName(); }
        static constexpr std::string_view fullyQualifiedName() { return minSINCE(); }
      };

      struct MAX_SINCE {
        typedef cond::Time_t type;
        static constexpr size_t size = 0;
        static constexpr std::string_view tableName() { return SINCE::tableName(); }
        static constexpr std::string_view fullyQualifiedName() { return maxSINCE(); }
      };

      class Table : public IIOVTable {
      public:
        explicit Table(coral::ISchema& schema);
        ~Table() override {}
        bool exists() override;
        void create() override;
        size_t getGroups(const std::string& tag,
                         const boost::posix_time::ptime& snapshotTime,
                         unsigned long long groupSize,
                         std::vector<cond::Time_t>& groups) override;
        size_t select(const std::string& tag,
                      cond::Time_t lowerGroup,
                      cond::Time_t upperGroup,
                      const boost::posix_time::ptime& snapshotTime,
                      std::vector<std::tuple<cond::Time_t, cond::Hash> >& iovs) override;
        bool getLastIov(const std::string& tag,
                        const boost::posix_time::ptime& snapshotTime,
                        cond::Time_t& since,
                        cond::Hash& hash) override;
        bool getSize(const std::string& tag, const boost::posix_time::ptime& snapshotTime, size_t& size) override;
        bool getRange(const std::string& tag,
                      cond::Time_t begin,
                      cond::Time_t end,
                      const boost::posix_time::ptime& snapshotTime,
                      std::vector<std::tuple<cond::Time_t, cond::Hash> >& iovs) override;
        void insertOne(const std::string& tag,
                       cond::Time_t since,
                       cond::Hash payloadHash,
                       const boost::posix_time::ptime& insertTime) override;
        void insertMany(
            const std::string& tag,
            const std::vector<std::tuple<cond::Time_t, cond::Hash, boost::posix_time::ptime> >& iovs) override;
        void eraseOne(const std::string& tag, cond::Time_t since, cond::Hash payloadId) override;
        void eraseMany(const std::string& tag, const std::vector<std::tuple<cond::Time_t, cond::Hash> >& iovs) override;
        void erase(const std::string& tag) override;

      private:
        coral::ISchema& m_schema;
      };
    }

    conddb_table(TAG_AUTHORIZATION) {
      conddb_column(TAG_NAME, std::string);
      conddb_column(ACCESS_TYPE, int);
      conddb_column(CREDENTIAL, std::string);
      conddb_column(CREDENTIAL_TYPE, int);

      class Table : public ITagAccessPermissionTable {
      public:
        explicit Table(coral::ISchema& schema);
        ~Table() override {}
        bool exists() override;
        void create() override;
        bool getAccessPermission(const std::string& tagName,
                                 const std::string& credential,
                                 int credentialType,
                                 int accessType) override;
        void setAccessPermission(const std::string& tagName,
                                 const std::string& credential,
                                 int credentialType,
                                 int accessType) override;
        void removeAccessPermission(const std::string& tagName,
                                    const std::string& credential,
                                    int credentialType) override;
        void removeEntriesForCredential(const std::string& credential, int credentialType) override;

      private:
        coral::ISchema& m_schema;
      };
    }

    conddb_table(TAG_LOG) {
      conddb_column(TAG_NAME, std::string);
      conddb_column(EVENT_TIME, boost::posix_time::ptime);
      conddb_column(USER_NAME, std::string);
      conddb_column(HOST_NAME, std::string);
      conddb_column(COMMAND, std::string);
      conddb_column(ACTION, std::string);
      conddb_column(USER_TEXT, std::string);

      class Table : public ITagLogTable {
      public:
        explicit Table(coral::ISchema& schema);
        ~Table() override {}
        bool exists() override;
        void create() override;
        void insert(const std::string& tag,
                    const boost::posix_time::ptime& eventTime,
                    const std::string& userName,
                    const std::string& hostName,
                    const std::string& command,
                    const std::string& action,
                    const std::string& userText) override;

      private:
        coral::ISchema& m_schema;
      };
    }

    class IOVSchema : public IIOVSchema {
    public:
      explicit IOVSchema(coral::ISchema& schema);
      ~IOVSchema() override {}
      bool exists() override;
      bool create() override;
      ITagTable& tagTable() override;
      IIOVTable& iovTable() override;
      ITagLogTable& tagLogTable() override;
      ITagAccessPermissionTable& tagAccessPermissionTable() override;
      IPayloadTable& payloadTable() override;

    private:
      TAG::Table m_tagTable;
      IOV::Table m_iovTable;
      TAG_LOG::Table m_tagLogTable;
      TAG_AUTHORIZATION::Table m_tagAccessPermissionTable;
      PAYLOAD::Table m_payloadTable;
    };

  }  // namespace persistency
}  // namespace cond
#endif
